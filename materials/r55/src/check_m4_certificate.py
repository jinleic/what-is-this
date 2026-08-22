#!/usr/bin/env python3
"""Independent checker for the m=4 higher-identity certificate artifact.

Role
----
`search_m4_cuts.py` writes `data/higher_identity_m4.json` (schema_version 2,
campaign higher_order_identity_positive_deficiency_m4) as the single
certificate artifact of the m=4 campaign.  This module re-verifies that
artifact from scratch, sharing no code with the producer chain: it imports
neither `subgraph_identities`, `m3_deficiency_cone`, `m4_subgraph_identities`,
`m4_deficiency_cone`, nor `search_m4_cuts`.  Only the frozen graph6 parser
`check_ramsey.parse_graph6_line` is imported; every kernel primitive,
selection rule, window, state, count, hash, certificate, witness, and
disposition below is recomputed here from the frozen constants and the
recorded input bytes, and any disagreement is a `CheckViolation` naming the
offending field.

The checker's own m=4 kernel (re-derived from scratch)
-----------------------------------------------------
* `t`: triangles by common-neighbour bit counts over every edge.
* `diamond`: for every non-edge `uv`, the adjacent pairs inside the common
  neighbourhood `A_u & A_v` -- each induced diamond is counted exactly once,
  by its unique non-edge (its degree-2 pair).
* `s(K3+K1)`: for every triangle `T`, the vertices adjacent to none of `T`
  (excluded `T` itself; the triangle's own vertices are killed by their
  mutual adjacency, so no explicit removal is needed).
* `i4`: `(1/6) * sum over non-edges uv of [C(c,2) - e(W)]` with `W` the
  common non-neighbourhood of `u` and `v` excluding both.  The edge
  subtraction inside `W` is REQUIRED: without it the sum equals `6*i4 + X`
  where `X` counts induced one-edge quartets (one-edge graph on six
  vertices: exact `i4 = 9`, shortcut `10`).
* `K4`: for every triangle, its common neighbours; every K4 is counted by
  its four triangles, so the accumulator divides by four.
* Universal residual (holds for EVERY graph; scaled by twelve):

      sum_v [ 3(n+2-2d_v) t(N_v) + 4 diamond(N_v) - 6 s(K3+K1,N_v)
              + 12 K4(N_v) - 12 K4(D_v) ] = 0

  where `N_v` is the neighbourhood and `D_v` the dual of `v`.  The
  `diamond_coef` and `k4_correction` keyword arguments exist so mutated
  coefficients can be shown to break the identity (K5 minus one edge has
  diamond neighbourhoods; K5 has K4 neighbourhoods).
* Specialized residual (graphs whose every neighbourhood is K4-free, in
  particular every (5,5)-graph): the same sum without the `+12 K4(N_v)`
  term; a K4 inside any neighbourhood raises `ValueError`.

What is checked, in fail-fast stage order
-----------------------------------------
1.  Schema preflight: exact key sets everywhere, canonical `r55/...`
    relative paths, no absolute paths anywhere in the document, sorted
    unique lists, per-section self-consistency (histogram sums, catalog
    windows against their own histograms, sound outer-window bounds,
    state invariants).
2.  Hash preflight: every recorded input exists under the root, its size
    equals the recorded byte count, and its streamed SHA-256 equals the
    recorded value.
3.  Validation selection and derived inputs: validation records are
    re-selected under the frozen rules (r44 decoys skipped, exactly 36
    extreme classes at orders 17..23, r35 window coverage 0..13), the
    recorded input set must equal the derived one, catalog line counts
    and recorded digests are tied back to their validation records.
4.  Catalog rescan: every extreme-class and census graph is re-parsed and
    re-counted with the checker's OWN five-motif kernel; the observed
    per-class histograms, owners, exact windows, and q outer intervals
    must equal the recorded ones, and every observed motif value must lie
    inside its recorded outer window.
5.  n=49 replay: the census is streamed, the 132-edge order-24 graphs are
    re-counted, and their t/diamond/s(K3+K1)/i4 aggregates must reproduce
    the frozen MR97 Thm 3.2 numbers (row 1584 = 12*132, forced mean 132,
    observed multiset [138, 144] strictly above the mean).
6.  State rebuild: windows, edge bounds, and the cone are rebuilt by the
    sole g formula and the sign-aware h formula; every recorded field,
    including the four provenance strings, must equal the rebuild.
7.  Search evidence: the search records are verified against the REBUILT
    states (not the stored ones), so a corrupted coefficient, bound,
    weight, or status is caught against independently derived truth.
8.  Replay: the 328 published R(5,5,42) graphs and their complements are
    streamed and the checker's OWN specialized m=4 residual kernel must
    return zero on every one; the replay section must record exactly what
    the file yields.
9.  Terminal mapping: the disposition must equal the verdict implied by
    the three route statuses (or `M4_ANALYZED` for the pre-search
    artifact with no searches).
10. Final re-hash: every recorded input is hashed again after the whole
    analysis, so a mid-run mutation of any file is still caught.

Outer-window policy (disclosed deviation): the producer's envelope windows
come from an exact optimization over the R1-R5 relation polytope whose
certified enumeration is budgeted at thirty minutes per run; recomputing
that LP inside every checker run would dominate the verification.  The
checker instead (a) recomputes the m=3 q outer intervals exactly from the
degree-histogram relaxation, (b) requires every recorded outer window to
stay inside the trivially sound combinatorial bounds, and (c) requires
every observed motif value of the class to lie inside the recorded window.
The h endpoints of missing strata are then resolved through these
recorded windows; present strata always use the re-streamed exact catalog
windows.  Scoping (review 2026-08-18): for PRESENT classes the audit is a
genuine soundness audit (caps plus observed-value containment); for
MISSING classes it is caps and schema only, so envelope soundness there
remains anchored to trust root 6 (the hash-pinned producer derivation).
The delivered M4_NO_CUT_IN_FROZEN_BASIS verdict is unaffected: its exact
witnesses prove non-provability over the recorded cone, and the g/balance
sides are recomputed exactly.  A future M4_ACCEPTED_CUT disposition is
NOT evidence-grade under this policy: it requires the plan-letter R1-R5
LP recompute per class before the verdict label holds.

CLI
---
`check_m4_certificate.py [--root DIR] ANALYSIS` exits 0 and prints
`M4 EVIDENCE VERIFIED: <disposition>` on the last stdout line, or exits 1
with `M4 EVIDENCE FAILED: <field>: <counterexample>` on stderr.  `--root`
substitutes the math repository root that every recorded `r55/...` path
resolves against; the default is the checker's own `parents[2]`.
"""

import argparse
import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import parse_graph6_line  # noqa: E402


# ------------------------------------------------------- frozen constants ---

N_TARGET = 45        # order of the hypothetical R(5,5,45) graph
R44 = 18             # R(4,4): lower bound of an R(4,5) degree window
R35_MAX_ORDER = 13   # R(3,5) - 1: upper bound of an R(4,5) degree window
R45_MAX_ORDER = 24   # R(4,5) - 1: maximum order of an R(4,5) graph
STATE_DEGREES = tuple(
    range(N_TARGET - 1 - R45_MAX_ORDER, R45_MAX_ORDER + 1)
)
VERTICES = Fraction(N_TARGET)   # every vertex occupies exactly one state

SCHEMA_VERSION = 2
CAMPAIGN_ID = "higher_order_identity_positive_deficiency_m4"

DATA_DIRNAME = "data"           # r55/data owns every catalog and table
TABLE_FILE = "structural_tables.json"
VALIDATION_FILE = "VALIDATION.json"
VALIDATION_EXTREME_FILE = "VALIDATION_extreme.json"
R35_FILES = tuple(
    f"r35_{order}.g6" for order in range(1, R35_MAX_ORDER + 1)
)
R45_CENSUS_FILE = "r45_24.g6"
CENSUS_ORDER = R45_MAX_ORDER
CENSUS_EDGE_BAND = (116, 132)   # published extremal edge band at order 24
R45_EXTREME_DIRNAME = "r45extreme"
R45_EXTREME_ORDERS = tuple(range(17, R45_MAX_ORDER))   # 17..23
R45_EXTREME_RECORDS = 36        # validated fixed-edge classes in that band
R55_REPLAY_FILE = "r55_42some.g6"
R55_REPLAY_ORDER = 42
# The executed sources pinned by the producing tasks; the artifact's input
# list must record exactly these four Python files (the frozen m=3
# machinery, the shared identity module, and the m=3 artifact itself may
# appear as additional hash-pinned provenance).
RECORDED_SOURCES = (
    "check_ramsey.py", "m4_deficiency_cone.py", "m4_subgraph_identities.py",
    "search_m4_cuts.py",
)
OPTIONAL_INPUTS = (
    "r55/src/m3_deficiency_cone.py",
    "r55/src/subgraph_identities.py",
    "r55/data/higher_identity_m3.json",
)

MOTIFS = ("q", "t", "diamond", "k3k1", "i4")
OUTER_SOURCE = "outer_degree_histogram"
TRIVIAL_MOTIF_SOURCE = "trivial:binomial_coefficient_bound"
ENVELOPE_SOURCE = "envelope_lp"
CATALOG_TAG = "catalog:"
OUTER_METHODS = (
    "catalog", "envelope_lp", "trivial_fallback", "trivial_sound_window",
)
TRUST_ROOT_COUNT = 6            # the plan's Trust Roots section

# The n=49 ground truth of McKay-Radziszowski Thm 3.2.
N49_ORDER = 49
N49_DEGREE = 24                 # a 24-regular hypothetical R(5,5,49) graph
N49_CENSUS_EDGES = 132          # extremal R(4,5,24) edge count

ACCEPTED = "ACCEPTED_EXACT_CUT"
REJECTED = "REJECTED_BY_EXACT_WITNESS"
UNRESOLVED = "CERTIFICATION_UNRESOLVED"
UNAVAILABLE = "UNAVAILABLE_NO_COVER_CERTIFICATE"
ROUTE_STATUSES = (ACCEPTED, REJECTED, UNRESOLVED)
ROUTE_IDS = ("total_deficiency", "degree20_count",
             "deficiency_ge8_count")
UNAVAILABLE_ROUTE = "required_local_family"

DISPOSITION_ACCEPTED = "M4_ACCEPTED_CUT"
DISPOSITION_NO_CUT = "M4_NO_CUT_IN_FROZEN_BASIS"
DISPOSITION_UNRESOLVED = "M4_CERTIFICATION_UNRESOLVED"
DISPOSITION_ANALYZED = "M4_ANALYZED"

UNAVAILABLE_RECORD = {
    "objective_id": UNAVAILABLE_ROUTE,
    "certificate": None,
    "exact_upper_bound": None,
    "primal_witness": None,
    "exact_witness_value": None,
    "route_status": UNAVAILABLE,
}

_DEFAULT_ROOT = Path(__file__).resolve().parents[2]
_STREAM_BLOCK = 1 << 20


class CheckViolation(ValueError):
    """One artifact field disagrees with the independently derived truth."""


class InfeasibleEdgeClass(ValueError):
    """No integer degree histogram survives for an (order, edges) class."""


# ------------------------------------------------------------- combinatorics


def _c2(n):
    return n * (n - 1) // 2


def _c3(n):
    return n * (n - 1) * (n - 2) // 6


def _c4(n):
    return n * (n - 1) * (n - 2) * (n - 3) // 24


def _trivial_cap(order, motif):
    """Sound upper bound of one motif aggregate over ALL order-r graphs.

    Both frozen window conventions must fit: the binomial window
    [0, C(r,4)] of the producer's fail-closed fallback and the fixture's
    q cap r*C(r-1,2); every other aggregate is bounded by C(r,4) (t by
    C(r,3) <= C(r,4) for r >= 7, and the relaxed polytope multiples of
    the R1-R5 system stay below it on every order this cone touches).
    """
    if motif == "q":
        return max(order * _c2(order - 1), _c4(order))
    return _c4(order)


# ------------------------------------------------------ kernel primitives ---
# Re-derived here from scratch; the producer chain is never imported.

def complement(adj):
    """Complement graph, same vertex order, no loops."""
    n = len(adj)
    full = (1 << n) - 1
    return [full & ~adj[v] & ~(1 << v) for v in range(n)]


def _pair_motifs(adj, mask):
    """(t, diamond, k3k1, k4) of the subgraph induced on the vertex mask.

    All rows stay in the global vertex numbering; neighbourhoods are
    handled by masking the adjacency rows, so no induced copy is built.
    Triangle bookkeeping: each triangle appears once (u < v < w), the k3k1
    counter reads the in-mask common non-neighbourhood of the triangle
    (its own vertices are excluded by their mutual adjacency), and the k4
    accumulator sees each K4 once per triangle (divided by four at the
    end).  Diamonds are counted per non-edge via the adjacent pairs
    inside the common neighbourhood (each pair seen twice, halved).
    """
    triangles = diamonds2 = k3k1 = k4acc = 0
    pending = mask
    while pending:
        low_u = pending & -pending
        u = low_u.bit_length() - 1
        pending ^= low_u
        row_u = adj[u] & mask
        # --- edges (u, v) with v > u inside the mask: triangles ---
        later = row_u & ~((1 << (u + 1)) - 1)
        while later:
            low_v = later & -later
            v = low_v.bit_length() - 1
            later ^= low_v
            row_v = adj[v] & mask
            common = row_u & row_v
            tri = common & ~((1 << (v + 1)) - 1)
            while tri:
                low_w = tri & -tri
                w = low_w.bit_length() - 1
                tri ^= low_w
                triangles += 1
                k3k1 += (
                    mask & ~row_u & ~row_v & ~adj[w]
                ).bit_count()
                k4acc += (common & adj[w]).bit_count()
        # --- non-edges (u, v) with v > u inside the mask: diamonds ---
        nonedge = mask & ~adj[u] & ~((1 << (u + 1)) - 1)
        while nonedge:
            low_v = nonedge & -nonedge
            v = low_v.bit_length() - 1
            nonedge ^= low_v
            share = row_u & (adj[v] & mask)
            inside = share
            while inside:
                low_w = inside & -inside
                w = low_w.bit_length() - 1
                inside ^= low_w
                diamonds2 += (share & adj[w]).bit_count()
    return triangles, diamonds2 // 2, k3k1, k4acc // 4


def _k4_of_mask(adj, mask):
    """K4 count of the subgraph induced on `mask`, triangles only."""
    acc = 0
    pending = mask
    while pending:
        low_u = pending & -pending
        u = low_u.bit_length() - 1
        pending ^= low_u
        row_u = adj[u] & mask
        later = row_u & ~((1 << (u + 1)) - 1)
        while later:
            low_v = later & -later
            v = low_v.bit_length() - 1
            later ^= low_v
            common = row_u & (adj[v] & mask)
            tri = common & ~((1 << (v + 1)) - 1)
            while tri:
                low_w = tri & -tri
                tri ^= low_w
                acc += (common & adj[low_w.bit_length() - 1]).bit_count()
    return acc // 4


def independent_quad_count(adj):
    """i4 = (1/6) sum over non-edges uv of [C(c,2) - e(W)], W the common
    non-neighbourhood of u and v excluding both.  The induced-edge
    subtraction inside W is what the research shortcut omits."""
    n = len(adj)
    full = (1 << n) - 1
    total = 0
    for u in range(n):
        row_u = adj[u]
        rest = full & ~row_u & ~((1 << (u + 1)) - 1)
        while rest:
            low_v = rest & -rest
            v = low_v.bit_length() - 1
            rest ^= low_v
            w = full & ~row_u & ~adj[v] & ~low_v & ~(1 << u)
            c = w.bit_count()
            total += _c2(c)
            edges2 = 0
            scan = w
            while scan:
                low_x = scan & -scan
                x = low_x.bit_length() - 1
                scan ^= low_x
                edges2 += (w & adj[x]).bit_count()
            total -= edges2 // 2
    if total % 6:
        raise CheckViolation(
            f"kernel: independent quad accumulator {total} is not a "
            "multiple of six"
        )
    return total // 6


def motif4_vector(adj):
    """The five streaming aggregates of one graph, from the own kernel."""
    n = len(adj)
    wedges = 0
    degree_sum = 0
    for row in adj:
        degree = row.bit_count()
        degree_sum += degree
        wedges += _c2(degree)
    if degree_sum % 2:
        raise CheckViolation("kernel: odd degree sum")
    triangles, diamonds, k3k1, _k4 = _pair_motifs(adj, (1 << n) - 1)
    return {
        "q": wedges - triangles,
        "t": triangles,
        "diamond": diamonds,
        "k3k1": k3k1,
        "i4": independent_quad_count(adj),
    }


def m4_residual_universal(adj, diamond_coef=4, k4_correction=12):
    """The checker's own universal m=4 residual, scaled by twelve.

    Sum over every vertex of 3(n+2-2d)t(N_v) + diamond_coef*diamond(N_v)
    - 6*s(K3+K1,N_v) + k4_correction*(K4(N_v) - K4(D_v)); exactly zero
    for EVERY graph under the frozen coefficients.  The two keyword
    coefficients exist so mutated identities can be shown to break.
    """
    n = len(adj)
    full = (1 << n) - 1
    total = 0
    for v in range(n):
        row = adj[v]
        dual = full & ~row & ~(1 << v)
        t, diamonds, k3k1, k4_x = _pair_motifs(adj, row)
        k4_y = _k4_of_mask(adj, dual)
        total += (
            3 * (n + 2 - 2 * row.bit_count()) * t
            + diamond_coef * diamonds
            - 6 * k3k1
            + k4_correction * (k4_x - k4_y)
        )
    return total


def m4_residual_specialized(adj):
    """The specialized residual; requires K4-free neighbourhoods.

    Raises ValueError when any neighbourhood carries a K4 (then the
    +12*K4(N_v) correction of the universal form cannot be dropped).
    """
    n = len(adj)
    full = (1 << n) - 1
    total = 0
    for v in range(n):
        row = adj[v]
        t, diamonds, k3k1, k4_x = _pair_motifs(adj, row)
        if k4_x:
            raise ValueError(
                "specialized identity needs K4-free neighborhoods"
            )
        k4_y = _k4_of_mask(adj, full & ~row & ~(1 << v))
        total += (
            3 * (n + 2 - 2 * row.bit_count()) * t
            + 4 * diamonds
            - 6 * k3k1
            - 12 * k4_y
        )
    return total


# --------------------------------------------------------- value guards ---

def _exact_int(value, label):
    """Reject anything that is not already an exact Python int."""
    if type(value) is not int:
        raise CheckViolation(f"{label}: {value!r} is not an exact integer")
    return value


def _canonical_fraction(text, label):
    """The value of an already reduced canonical fraction string.

    One round trip rejects padded or signed text, decimal strings,
    unreduced fractions, and `n/1` spellings in a single stroke.
    """
    if type(text) is not str:
        raise CheckViolation(f"{label}: {text!r} must be a string")
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError):
        raise CheckViolation(
            f"{label}: {text!r} is not a rational"
        ) from None
    if str(value) != text:
        raise CheckViolation(
            f"{label}: {text!r} is not canonical, want {value}"
        )
    return value


def _require_keys(record, keys, label):
    """Exactly the schema keys, nothing missing and nothing extra."""
    if not isinstance(record, dict):
        raise CheckViolation(f"{label}: not an object")
    if tuple(sorted(record)) != tuple(keys):
        raise CheckViolation(
            f"{label}: has keys {sorted(record)}, expected {list(keys)}"
        )
    return record


def _is_sha256(value):
    """True only for a lowercase 64-hex digest string."""
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(char in "0123456789abcdef" for char in value)
    )


def _iter_strings(value):
    """Every string in the document, keys included."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield key
            for text in _iter_strings(item):
                yield text
    elif isinstance(value, list):
        for item in value:
            for text in _iter_strings(item):
                yield text


def _load_json(path, label):
    """One parsed JSON object, or a field-naming violation."""
    try:
        raw = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise CheckViolation(
            f"{label}: unreadable JSON ({exc})"
        ) from None
    if not isinstance(raw, dict):
        raise CheckViolation(f"{label}: top-level JSON object expected")
    return raw


def _sha256_of(path):
    """Streamed SHA-256 of one recorded input file."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(_STREAM_BLOCK), b""):
                digest.update(block)
    except OSError as exc:
        raise CheckViolation(
            f"sha256: {path.name}: unreadable input ({exc})"
        ) from None
    return digest.hexdigest()


# ---------------------------------------------------- validation selection ---

def _record_filename(record, source):
    """Filename and parent directory name of a validation record.

    The recorded `file` path is an absolute machine path of the producing
    run, so only its last two components are ever consumed.
    """
    raw = record.get("file")
    if not isinstance(raw, str) or not raw:
        raise CheckViolation(f"{source}: result record without a file field")
    pure = PurePosixPath(raw)
    if pure.name in ("", ".", ".."):
        raise CheckViolation(f"{source}: record file {raw!r} has no filename")
    return pure.name, pure.parent.name


def _require_validated(record, keys, label, source):
    """Every named predicate exactly true, and no recorded violation."""
    for key in keys:
        if record.get(key) is not True:
            raise CheckViolation(
                f"{source}: {label} does not have {key} == true"
            )
    violations = record.get("violations")
    if not isinstance(violations, list) or violations:
        raise CheckViolation(
            f"{source}: {label} has violations {violations!r}"
        )


def _record_order_and_count(record, label, source):
    """The single order of a record and its agreed graph count."""
    orders = record.get("n_values")
    if not isinstance(orders, dict) or len(orders) != 1:
        raise CheckViolation(
            f"{source}: {label} does not cover exactly one order"
        )
    key, count = next(iter(orders.items()))
    if not (isinstance(key, str) and key.isascii() and key.isdigit()):
        raise CheckViolation(
            f"{source}: {label} has a non-numeric order {key!r}"
        )
    order = int(key)
    count = _exact_int(count, f"{label} n_values[{order}]")
    graphs = _exact_int(record.get("graphs"), f"{label} graphs")
    if graphs <= 0 or graphs != count:
        raise CheckViolation(
            f"{source}: {label} claims {graphs} graphs but {count} at "
            f"order {order}"
        )
    return order, graphs


def select_main_records(path):
    """The validated records of the catalogs, census, and replay file.

    Wanted files: the 13 R(3,5) catalogs, the order-24 census, and the
    R(5,5,42) replay file, all under `r55/data/`.  r44_* records in the
    same report are (4,4) decoys and are skipped.  Graph counts are read
    from the report, never pinned to published values: the rescan must
    reproduce them from the bytes.
    """
    source = path.name
    raw = _load_json(path, source)
    results = raw.get("results")
    if not isinstance(results, list):
        raise CheckViolation(f"{source}: no results list")
    wanted = {
        name: (3, 5, order) for order, name in enumerate(R35_FILES, start=1)
    }
    wanted[R45_CENSUS_FILE] = (4, 5, CENSUS_ORDER)
    wanted[R55_REPLAY_FILE] = (5, 5, R55_REPLAY_ORDER)
    selected = {}
    for record in results:
        if not isinstance(record, dict):
            raise CheckViolation(f"{source}: malformed result record")
        name, parent = _record_filename(record, source)
        if name not in wanted:
            continue            # the (4,4) decoys and anything unknown
        if parent != DATA_DIRNAME:
            raise CheckViolation(
                f"{source}: {name} is recorded outside {DATA_DIRNAME}/"
            )
        if name in selected:
            raise CheckViolation(
                f"duplicate: {source}: duplicate record for {name}"
            )
        want_s, want_t, want_order = wanted[name]
        if (record.get("s"), record.get("t")) != (want_s, want_t):
            raise CheckViolation(
                f"{source}: {name} is not an R({want_s},{want_t}) record"
            )
        _require_validated(
            record, ("ok", "count_matches", "all_ramsey"), name, source
        )
        order, graphs = _record_order_and_count(record, name, source)
        if order != want_order:
            raise CheckViolation(
                f"{source}: {name} covers order {order}, expected "
                f"{want_order}"
            )
        digest = record.get("sha256")
        if not _is_sha256(digest):
            raise CheckViolation(f"{source}: {name} has no recorded sha256")
        selected[name] = {"order": order, "graphs": graphs, "sha256": digest}
    missing = sorted(set(wanted) - set(selected))
    if missing:
        raise CheckViolation(
            f"{source}: no validated record for {missing}"
        )
    return selected


def select_r35_windows(path):
    """Complete R(3,5) edge windows, order zero synthesized as (0, 0)."""
    source = path.name
    raw = _load_json(path, source)
    results = raw.get("results")
    if not isinstance(results, list):
        raise CheckViolation(f"{source}: no results list")
    windows = {0: (0, 0)}       # the empty graph: no vertices, no edges
    for record in results:
        if not isinstance(record, dict):
            raise CheckViolation(f"{source}: malformed result record")
        if (record.get("s"), record.get("t")) != (3, 5):
            continue
        if not (record.get("ok") and record.get("count_matches")):
            continue            # unvalidated records cannot define a window
        order, _ = _record_order_and_count(record, "R(3,5) record", source)
        low = _exact_int(record.get("edge_min"), f"R(3,5,{order}) edge_min")
        high = _exact_int(record.get("edge_max"), f"R(3,5,{order}) edge_max")
        if low > high:
            raise CheckViolation(
                f"r35_edge_windows: order {order} window is empty "
                f"({low} > {high})"
            )
        known = windows.setdefault(order, (low, high))
        if known != (low, high):
            raise CheckViolation(
                f"r35_edge_windows: conflicting windows for order {order}: "
                f"{known} and {(low, high)}"
            )
    missing = [
        order for order in range(1, R35_MAX_ORDER + 1) if order not in windows
    ]
    if missing:
        raise CheckViolation(
            f"r35_edge_windows: no validated R(3,5) window for orders "
            f"{missing}"
        )
    return windows


def _extreme_name_class(name, source):
    """The (order, edges) class encoded by an `r45<order>.<edges>.g6` name."""
    if not name.startswith("r45") or not name.endswith(".g6"):
        raise CheckViolation(
            f"{source}: {name} is not an r45 edge-class catalog"
        )
    fields = name[3:-3].split(".")
    if len(fields) != 2 or not all(
        field.isascii() and field.isdigit() for field in fields
    ):
        raise CheckViolation(
            f"{source}: {name} does not encode order and edges"
        )
    return int(fields[0]), int(fields[1])


def select_extreme_records(path):
    """Every validated fixed-edge R(4,5) class at orders 17..23.

    Exactly 36 classes, every order of the band covered.  This report
    must carry no `count_matches` or `sha256` field: absence is
    required, not tolerated; filename, metadata order, and metadata
    edge class must agree instead.
    """
    source = path.name
    raw = _load_json(path, source)
    results = raw.get("results")
    if not isinstance(results, list):
        raise CheckViolation(f"{source}: no results list")
    selected = {}
    for record in results:
        if not isinstance(record, dict):
            raise CheckViolation(f"{source}: malformed result record")
        name, parent = _record_filename(record, source)
        order, graphs = _record_order_and_count(record, name, source)
        if order not in R45_EXTREME_ORDERS:
            continue            # the archive holds orders 4..16 as well
        if parent != R45_EXTREME_DIRNAME:
            raise CheckViolation(
                f"{source}: target-order record {name} is outside "
                f"{R45_EXTREME_DIRNAME}/"
            )
        if (record.get("s"), record.get("t")) != (4, 5):
            raise CheckViolation(f"{source}: {name} is not an R(4,5) record")
        _require_validated(record, ("ok", "all_ramsey"), name, source)
        for field in ("count_matches", "sha256"):
            if field in record:
                raise CheckViolation(
                    f"{source}: {name} carries a forbidden {field} field"
                )
        edge_min = _exact_int(record.get("edge_min"), f"{name} edge_min")
        edge_max = _exact_int(record.get("edge_max"), f"{name} edge_max")
        if edge_min != edge_max:
            raise CheckViolation(
                f"{source}: {name} spans edges {edge_min}..{edge_max}, not "
                "one fixed class"
            )
        if _extreme_name_class(name, source) != (order, edge_min):
            raise CheckViolation(
                f"{source}: {name} disagrees with its metadata "
                f"(order {order}, edges {edge_min})"
            )
        if name in selected:
            raise CheckViolation(
                f"duplicate: {source}: duplicate record for {name}"
            )
        selected[name] = {"order": order, "edges": edge_min, "graphs": graphs}
    if len(selected) != R45_EXTREME_RECORDS:
        raise CheckViolation(
            f"extreme: {len(selected)} validated order-"
            f"{R45_EXTREME_ORDERS[0]}..{R45_EXTREME_ORDERS[-1]} records in "
            f"{source}, expected {R45_EXTREME_RECORDS}"
        )
    covered = {info["order"] for info in selected.values()}
    missing = [order for order in R45_EXTREME_ORDERS if order not in covered]
    if missing:
        raise CheckViolation(
            f"extreme: no validated edge class for orders {missing}"
        )
    return selected


def load_edge_bounds(path):
    """Published R(4,5) edge extrema per order, from the owning table."""
    raw = _load_json(path, TABLE_FILE)
    table = raw.get("extremal_edge_tables_R45")
    if not isinstance(table, dict):
        raise CheckViolation(
            f"{TABLE_FILE}: no extremal_edge_tables_R45 object"
        )
    bounds = {}
    for order, values in table.items():
        if not isinstance(values, dict):
            raise CheckViolation(
                f"{TABLE_FILE}: order {order!r} is not an object"
            )
        low = _exact_int(values.get("min"), f"order {order} min")
        high = _exact_int(values.get("max"), f"order {order} max")
        if low > high:
            raise CheckViolation(
                f"{TABLE_FILE}: order {order} window is empty ({low} > {high})"
            )
        bounds[int(order)] = (low, high)
    missing = [order for order in STATE_DEGREES if order not in bounds]
    dual_missing = [
        order for order in (N_TARGET - 1 - d for d in STATE_DEGREES)
        if order not in bounds
    ]
    missing = sorted(set(missing) | set(dual_missing))
    if missing:
        raise CheckViolation(
            f"{TABLE_FILE}: no edge extrema for cone orders {missing}"
        )
    return bounds


# ------------------------------------------------- envelope and cone math ---

def degree_bounds(order):
    """Ramsey degree window of an R(4,5,order) graph."""
    return max(0, order - R44), min(R35_MAX_ORDER, order - 1)


def _iter_q_terms(low, high, vertices, degree_sum, wedges, lows, highs):
    """Every degree histogram's (wedge, R(3,5) lo, hi) aggregates.

    Counts range over degrees low..high with fixed vertex and degree
    sums; graphicality is deliberately not imposed: the family is an
    outer relaxation, so skipping it only widens the envelope.  The
    admissible counts at a level form one closed-form interval (with
    ascending degrees the floor condition excludes small counts and the
    cap and nonnegativity conditions exclude large ones), the three
    aggregates are carried incrementally, and the walk is an explicit
    stack instead of a recursive generator chain.  The yielded leaf set
    is exactly the per-histogram dot products of the frozen m=3
    enumeration.
    """
    if high < low or vertices < 0 or degree_sum < 0:
        return
    degrees = tuple(range(low, high + 1))
    top = degrees[-1]
    last = len(degrees) - 1
    if last == 0:
        if degrees[0] * vertices == degree_sum:
            yield (vertices * wedges[0], vertices * lows[0],
                   vertices * highs[0])
        return

    def count_range(degree, floor_degree, verts_left, sum_left):
        count_min = max(
            0,
            -((sum_left - verts_left * floor_degree)
              // (floor_degree - degree)),
        )
        if degree:
            count_max = min(
                verts_left,
                sum_left // degree,
                (verts_left * top - sum_left) // (top - degree),
            )
        else:
            count_max = min(
                verts_left,
                (verts_left * top - sum_left) // top,
            )
        return count_min, count_max

    frame_verts = [0] * last
    frame_sums = [0] * last
    frame_wedge = [0] * last
    frame_lo = [0] * last
    frame_hi = [0] * last
    frame_count = [0] * last
    frame_max = [0] * last

    count_min, count_max = count_range(
        degrees[0], degrees[1], vertices, degree_sum
    )
    if count_min > count_max:
        return
    frame_verts[0] = vertices
    frame_sums[0] = degree_sum
    frame_wedge[0] = frame_lo[0] = frame_hi[0] = 0
    frame_count[0] = count_min - 1
    frame_max[0] = count_max

    position = 0
    leaf_wedge = wedges[last]
    leaf_low = lows[last]
    leaf_high = highs[last]
    leaf_degree = degrees[last]
    while position >= 0:
        count = frame_count[position] + 1
        if count > frame_max[position]:
            position -= 1
            continue
        frame_count[position] = count
        degree = degrees[position]
        verts_rest = frame_verts[position] - count
        sum_rest = frame_sums[position] - count * degree
        wedge = frame_wedge[position] + count * wedges[position]
        lo_acc = frame_lo[position] + count * lows[position]
        hi_acc = frame_hi[position] + count * highs[position]
        if position == last - 1:
            if leaf_degree * verts_rest == sum_rest:
                yield (wedge + verts_rest * leaf_wedge,
                       lo_acc + verts_rest * leaf_low,
                       hi_acc + verts_rest * leaf_high)
            continue
        count_min, count_max = count_range(
            degrees[position + 1], degrees[position + 2],
            verts_rest, sum_rest,
        )
        if count_min > count_max:
            continue
        position += 1
        frame_verts[position] = verts_rest
        frame_sums[position] = sum_rest
        frame_wedge[position] = wedge
        frame_lo[position] = lo_acc
        frame_hi[position] = hi_acc
        frame_count[position] = count_min - 1
        frame_max[position] = count_max


def q_outer_interval(order, edges, windows):
    """Outer interval for q = sum_x binom(d_x, 2) - t over an edge class.

    Every integer degree histogram of the class is enumerated; for each,
    the R(3,5) windows bound 3t to [T3lo, T3hi], integrality of t refines
    that to [ceil(T3lo/3), floor(T3hi/3)], and the local q interval is
    [W - t_hi, W - t_lo].  The envelope is the union over histograms.
    """
    if order < 0 or edges < 0:
        raise InfeasibleEdgeClass(
            f"negative class (order={order}, edges={edges})"
        )
    low, high = degree_bounds(order)
    degrees = range(low, high + 1)
    missing = [degree for degree in degrees if degree not in windows]
    if missing:
        raise CheckViolation(
            f"r35_edge_windows: no R(3,5) window for degrees {missing} "
            f"(order={order}, edges={edges})"
        )
    wedges = [degree * (degree - 1) // 2 for degree in degrees]
    tri3_lo = [windows[degree][0] for degree in degrees]
    tri3_hi = [windows[degree][1] for degree in degrees]
    q_lo = None
    q_hi = None
    standing = 0
    for wedge_sum, lo_sum, hi_sum in _iter_q_terms(
        low, high, order, 2 * edges, wedges, tri3_lo, tri3_hi
    ):
        standing += 1
        t_lo = -(-lo_sum // 3)  # ceil towards the sound interior
        t_hi = hi_sum // 3
        if t_lo > t_hi:
            continue
        local_lo = wedge_sum - t_hi
        local_hi = wedge_sum - t_lo
        q_lo = local_lo if q_lo is None else min(q_lo, local_lo)
        q_hi = local_hi if q_hi is None else max(q_hi, local_hi)
    if q_lo is None:
        if standing:
            raise InfeasibleEdgeClass(
                f"no histogram meets the R(3,5) triangle window for "
                f"order={order}, edges={edges}"
            )
        raise InfeasibleEdgeClass(
            f"no degree histogram for order={order}, edges={edges}"
        )
    return q_lo, q_hi


def budget2(d, edge_bounds):
    """Twice the largest m=2 excess deficit a degree-d vertex can carry."""
    m = N_TARGET - 1 - d
    h2_min = (
        2 * (m * (m - 1) // 2 - edge_bounds[m][1] - edge_bounds[d][1])
        - d * (N_TARGET - 2 * d)
    )
    return -h2_min


def _cone_state(d, a, b, ex, ey, qx, qy, budget):
    """One cone cell; the sole g formula with the frozen sign pairing."""
    m = N_TARGET - 1 - d
    dual_base = 3 * (m * (m - 1) * (m - 2) // 6 - ey * (m - 2))
    rhs_base = (N_TARGET + 3 - 3 * d) * ex
    return {
        "d": d,
        "a": a,
        "b": b,
        "deficiency": a + b,
        "excess_balance": 2 * (a + b) - budget,
        "g_lo": dual_base + 3 * qy[0] - rhs_base - 3 * qx[1],
        "g_hi": dual_base + 3 * qy[1] - rhs_base - 3 * qx[0],
    }


def _h_endpoints(d, x_windows, y_windows):
    """The sign-aware m=4 row interval over the two strata windows.

    row = 3(47-2d)*t(x) + 4*diamond(x) - 6*k3k1(x) - 12*i4(y); the
    t coefficient 3(47-2d) flips sign at d = 24, so the t window is
    scaled min/max-aware while the other three coefficients are fixed
    sign.
    """
    coef = 3 * (N_TARGET + 2 - 2 * d)
    t_lo, t_hi = x_windows["t"]
    if coef >= 0:
        lo_t, hi_t = coef * t_lo, coef * t_hi
    else:
        lo_t, hi_t = coef * t_hi, coef * t_lo
    return (
        lo_t + 4 * x_windows["diamond"][0]
        - 6 * x_windows["k3k1"][1] - 12 * y_windows["i4"][1],
        hi_t + 4 * x_windows["diamond"][1]
        - 6 * x_windows["k3k1"][0] - 12 * y_windows["i4"][0],
    )


_MISSING_MOTIF_SOURCES = frozenset((TRIVIAL_MOTIF_SOURCE, ENVELOPE_SOURCE))


def rebuild_states(windows, edge_bounds, catalog_windows, class_source,
                   motif_windows, outer_windows):
    """Every state of the n=45 cone, sorted by (d, a, b).

    The g side follows the frozen m=3 rule: the exact catalog q window
    when the class is owned, else the degree-histogram relaxation; a
    pair is omitted only when a side class has no surviving histogram.
    The h side resolves exact catalog motif windows for owned classes
    and the soundness-audited recorded outer windows otherwise.
    """
    cache = {}

    def resolve(order, edges):
        key = (order, edges)
        if key not in cache:
            if key in catalog_windows:
                cache[key] = catalog_windows[key]
            else:
                try:
                    cache[key] = q_outer_interval(order, edges, windows)
                except InfeasibleEdgeClass:
                    cache[key] = None
        return cache[key]

    def motif_window(cls):
        if cls in motif_windows:
            return motif_windows[cls]
        if cls in outer_windows:
            return outer_windows[cls]
        raise CheckViolation(
            f"outer_motif_windows: no certified window for the stratum "
            f"class (order={cls[0]}, edges={cls[1]})"
        )

    states = []
    for d in STATE_DEGREES:
        m = N_TARGET - 1 - d
        d_low, d_high = edge_bounds[d]
        m_low, m_high = edge_bounds[m]
        budget = budget2(d, edge_bounds)
        for a in range(d_high - d_low + 1):
            ex = d_high - a
            qx = resolve(d, ex)
            if qx is None:
                continue
            x_windows = motif_window((d, ex))
            x_owner = class_source.get((d, ex))
            for b in range(m_high - m_low + 1):
                ey = m_high - b
                qy = resolve(m, ey)
                if qy is None:
                    continue
                y_windows = motif_window((m, ey))
                y_owner = class_source.get((m, ey))
                state = _cone_state(d, a, b, ex, ey, qx, qy, budget)
                state["h_lo"], state["h_hi"] = _h_endpoints(
                    d, x_windows, y_windows
                )
                state["x_interval_source"] = (
                    x_owner if x_owner else OUTER_SOURCE
                )
                state["y_interval_source"] = (
                    y_owner if y_owner else OUTER_SOURCE
                )
                state["x_motif_ok"] = (
                    frozenset((x_owner, CATALOG_TAG + x_owner))
                    if x_owner else _MISSING_MOTIF_SOURCES
                )
                state["y_motif_ok"] = (
                    frozenset((y_owner, CATALOG_TAG + y_owner))
                    if y_owner else _MISSING_MOTIF_SOURCES
                )
                states.append(state)
    states.sort(key=lambda state: (state["d"], state["a"], state["b"]))
    return states


# ---------------------------------------------------------- catalog rescan ---

def _stream_motif_catalog(path, expected_order, expected_edges,
                          expected_count):
    """One pass over a graph6 catalog: the five motif frequencies.

    Every non-empty line is parsed, its order checked, its edge count
    checked when the class is fixed, and the five aggregates of the
    checker's own kernel accumulated per (order, edges).  No graph is
    retained; memory stays proportional to the tiny frequency maps.
    """
    label = path.name
    counts = {}
    graphs = 0
    try:
        stream = path.open("r", encoding="ascii")
    except OSError as exc:
        raise CheckViolation(
            f"catalog: {label}: unreadable catalog ({exc})"
        ) from None
    with stream:
        try:
            for number, line in enumerate(stream, start=1):
                text = line.strip()
                if not text:
                    continue
                graphs += 1
                try:
                    order, adj = parse_graph6_line(text)
                except ValueError as exc:
                    raise CheckViolation(
                        f"catalog: {label} line {number}: bad graph6 ({exc})"
                    ) from None
                if order != expected_order:
                    raise CheckViolation(
                        f"catalog: {label} line {number}: order {order}, "
                        f"expected {expected_order}"
                    )
                edges = sum(row.bit_count() for row in adj) // 2
                if expected_edges is not None and edges != expected_edges:
                    raise CheckViolation(
                        f"catalog: {label} line {number}: {edges} edges, "
                        f"expected {expected_edges}"
                    )
                vector = motif4_vector(adj)
                frequency = counts.get((order, edges))
                if frequency is None:
                    frequency = counts[(order, edges)] = {}
                for motif, value in vector.items():
                    tally = frequency.get(motif)
                    if tally is None:
                        tally = frequency[motif] = {}
                    tally[value] = tally.get(value, 0) + 1
        except (OSError, UnicodeError) as exc:
            raise CheckViolation(
                f"catalog: {label}: unreadable catalog ({exc})"
            ) from None
    if not graphs:
        raise CheckViolation(f"catalog: {label}: empty catalog")
    if graphs != expected_count:
        raise CheckViolation(
            f"count: {label} holds {graphs} graphs, the validation report "
            f"records {expected_count}"
        )
    return graphs, counts


def sweep_catalogs(records):
    """Stream every catalog once: class histograms, line counts, owners.

    `records` is an iterable of `(path, expected_order,
    expected_edges_or_none, expected_graph_count, relative_path)`.  A
    class produced by two different files is an ownership violation.
    """
    histograms = {}
    owners = {}
    graph_counts = {}
    for path, order, edges, count, relative in records:
        location = path.as_posix()
        if location in graph_counts:
            raise CheckViolation(
                f"catalog: {path.name} is listed twice"
            )
        graphs, produced = _stream_motif_catalog(path, order, edges, count)
        graph_counts[location] = graphs
        for key, frequency in produced.items():
            owner = owners.get(key)
            if owner is not None:
                raise CheckViolation(
                    f"owner: class (order={key[0]}, edges={key[1]}) is "
                    f"produced by both {owner} and {relative}"
                )
            owners[key] = relative
            histograms[key] = {
                motif: {value: frequency[motif][value]
                        for value in sorted(frequency[motif])}
                for motif in MOTIFS
            }
    return histograms, owners, graph_counts


def count_graph6_lines(path):
    """Non-empty graph6 line count of one small catalog."""
    count = 0
    try:
        with path.open("r", encoding="ascii") as stream:
            for line in stream:
                if line.strip():
                    count += 1
    except (OSError, UnicodeError) as exc:
        raise CheckViolation(
            f"count: {path.name}: unreadable catalog ({exc})"
        ) from None
    if not count:
        raise CheckViolation(f"count: {path.name}: empty catalog")
    return count


# --------------------------------------------------------------- n49 replay ---

def _check_n49(document, root):
    """Re-derive the MR97 Thm 3.2 anchor from the census's own bytes."""
    n49 = document["n49_replay"]
    path = root / "r55" / DATA_DIRNAME / R45_CENSUS_FILE
    name = path.name
    row_coef = 3 * (N49_ORDER + 2 - 2 * N49_DEGREE)
    extremal = []
    try:
        stream = path.open("r", encoding="ascii")
    except OSError as exc:
        raise CheckViolation(
            f"n49_replay: {name}: unreadable census ({exc})"
        ) from None
    with stream:
        try:
            for number, line in enumerate(stream, start=1):
                text = line.strip()
                if not text:
                    continue
                try:
                    order, adj = parse_graph6_line(text)
                except ValueError as exc:
                    raise CheckViolation(
                        f"n49_replay: {name} line {number}: bad graph6 "
                        f"({exc})"
                    ) from None
                if order != CENSUS_ORDER:
                    raise CheckViolation(
                        f"n49_replay: {name} line {number}: order {order}, "
                        f"expected {CENSUS_ORDER}"
                    )
                edges = sum(row.bit_count() for row in adj) // 2
                if edges == N49_CENSUS_EDGES:
                    vector = motif4_vector(adj)
                    extremal.append((
                        row_coef * vector["t"]
                        + 4 * vector["diamond"]
                        - 6 * vector["k3k1"],
                        vector["i4"],
                    ))
        except (OSError, UnicodeError) as exc:
            raise CheckViolation(
                f"n49_replay: {name}: unreadable census ({exc})"
            ) from None
    recorded_count = n49["extremal_graphs"]
    if recorded_count != len(extremal):
        raise CheckViolation(
            f"extremal_graphs: {recorded_count} recorded, the census holds "
            f"{len(extremal)} {N49_CENSUS_EDGES}-edge order-"
            f"{CENSUS_ORDER} graphs"
        )
    rows = [row for row, _i4 in extremal]
    if any(row != rows[0] for row in rows):
        raise CheckViolation(
            f"row_scaled: the extremal graphs carry different rows {rows}"
        )
    if rows[0] != n49["row_scaled"]:
        raise CheckViolation(
            f"row_scaled: {n49['row_scaled']} recorded, the extremal "
            f"graphs give {rows[0]}"
        )
    if rows[0] % 12:
        raise CheckViolation(
            f"forced_i4_mean: row {rows[0]} is not a multiple of twelve"
        )
    mean = rows[0] // 12
    if mean != n49["forced_i4_mean"]:
        raise CheckViolation(
            f"forced_i4_mean: {n49['forced_i4_mean']} recorded, the row "
            f"forces {mean}"
        )
    observed = sorted(i4 for _row, i4 in extremal)
    if list(n49["observed_i4"]) != observed:
        raise CheckViolation(
            f"observed_i4: {n49['observed_i4']!r} recorded, the extremal "
            f"graphs give {observed}"
        )
    if observed[0] <= mean:
        raise CheckViolation(
            f"observed_i4: the smallest admissible {observed[0]} does not "
            f"exceed the forced mean {mean}, the n=49 refutation is gone"
        )


# ------------------------------------------------------------------ replay ---

def _check_replay(document, root):
    """Re-derive the replay section from the published file's own bytes.

    Every published graph and its complement must satisfy the checker's
    specialized residual kernel exactly, and the recorded counts must be
    exactly what the file yields: lines, lines, and twice lines.
    """
    replay = document["replay"]
    path = root / "r55" / DATA_DIRNAME / R55_REPLAY_FILE
    name = path.name
    lines = 0
    zeros = 0
    try:
        stream = path.open("r", encoding="ascii")
    except OSError as exc:
        raise CheckViolation(
            f"replay: {name}: unreadable replay file ({exc})"
        ) from None
    with stream:
        try:
            for number, line in enumerate(stream, start=1):
                text = line.strip()
                if not text:
                    continue
                lines += 1
                try:
                    order, adj = parse_graph6_line(text)
                except ValueError as exc:
                    raise CheckViolation(
                        f"replay: {name} line {number}: bad graph6 ({exc})"
                    ) from None
                if order != R55_REPLAY_ORDER:
                    raise CheckViolation(
                        f"replay: {name} line {number}: order {order}, "
                        f"expected {R55_REPLAY_ORDER}"
                    )
                for label, graph in (
                    ("published graph", adj),
                    ("complement", complement(adj)),
                ):
                    try:
                        residual = m4_residual_specialized(graph)
                    except ValueError as exc:
                        raise CheckViolation(
                            f"residual: {name} line {number}: the {label} "
                            f"breaks the specialized kernel ({exc})"
                        ) from None
                    if residual != 0:
                        raise CheckViolation(
                            f"residual: {name} line {number}: the {label} "
                            f"has residual {residual}, expected 0"
                        )
                    zeros += 1
        except (OSError, UnicodeError) as exc:
            raise CheckViolation(
                f"replay: {name}: unreadable replay file ({exc})"
            ) from None
    if replay["published_graphs"] != lines:
        raise CheckViolation(
            f"replay: published_graphs is {replay['published_graphs']}, "
            f"the file holds {lines}"
        )
    if replay["complements"] != lines:
        raise CheckViolation(
            f"replay: complements is {replay['complements']}, the file "
            f"holds {lines}"
        )
    if replay["residual_zero"] != 2 * lines:
        raise CheckViolation(
            f"replay residual_zero: {replay['residual_zero']} recorded, "
            f"the replay verifies {2 * lines}"
        )
    return lines, zeros


# -------------------------------------------------------- search evidence ---

_ROUTES = {
    # objective_id: (per-state value, frozen acceptance edge, inclusive)
    "total_deficiency": (
        lambda state: state["deficiency"], Fraction(315), True,
    ),
    "degree20_count": (
        lambda state: int(state["d"] == 20), Fraction(1), False,
    ),
    "deficiency_ge8_count": (
        lambda state: int(state["deficiency"] >= 8), Fraction(1), False,
    ),
}


def _route(objective_id):
    """The frozen route of an objective; an unregistered ID is rejected."""
    route = _ROUTES.get(objective_id)
    if route is None:
        raise CheckViolation(
            f"objective_id: {objective_id!r} is not a frozen route"
        )
    return route


def _clears_edge(value, limit, inclusive):
    """True when an exact value clears the frozen acceptance edge."""
    return value <= limit if inclusive else value < limit


def _verify_certificate(states, entry, objective_id):
    """Check one six-coefficient affine certificate exactly.

    Requires canonical coefficients with gamma, delta, epsilon, zeta >= 0
    and the inequality alpha + beta*xb + gamma*g_lo - delta*g_hi
    + epsilon*h_lo - zeta*h_hi - objective >= 0 at every state.
    Nothing here trusts how the coefficients were found.  Returns the
    exact 45*alpha.
    """
    _require_keys(entry, _CERTIFICATE_KEYS, "certificate")
    if entry["objective_id"] != objective_id:
        raise CheckViolation(
            f"certificate: objective {entry['objective_id']!r} is not the "
            f"record objective {objective_id!r}"
        )
    alpha = _canonical_fraction(entry["alpha"], "certificate alpha")
    beta = _canonical_fraction(entry["beta"], "certificate beta")
    gamma = _canonical_fraction(entry["gamma"], "certificate gamma")
    delta = _canonical_fraction(entry["delta"], "certificate delta")
    epsilon = _canonical_fraction(entry["epsilon"], "certificate epsilon")
    zeta = _canonical_fraction(entry["zeta"], "certificate zeta")
    for name, value in (("gamma", gamma), ("delta", delta),
                        ("epsilon", epsilon), ("zeta", zeta)):
        if value < 0:
            raise CheckViolation(
                f"certificate {name}: {value} is negative"
            )
    value = _route(objective_id)[0]
    for index, state in enumerate(states):
        slack = (
            alpha
            + beta * state["excess_balance"]
            + gamma * state["g_lo"]
            - delta * state["g_hi"]
            + epsilon * state["h_lo"]
            - zeta * state["h_hi"]
            - value(state)
        )
        if slack < 0:
            raise CheckViolation(
                f"certificate: state {index} (d={state['d']}, "
                f"a={state['a']}, b={state['b']}) violates by {-slack} "
                f"under alpha {alpha}, beta {beta}, gamma {gamma}, "
                f"delta {delta}, epsilon {epsilon}, zeta {zeta}"
            )
    return VERTICES * alpha


def _verify_witness(states, entry, objective_id):
    """Check one primal witness exactly; return its exact objective value.

    Strictly increasing in-range indices, canonical nonnegative weights,
    total weight 45, zero aggregate excess balance, and the three
    straddle pairs aggregate g_lo <= 0 <= g_hi and h_lo <= 0 <= h_hi:
    exactly the relations a real 45-vertex graph satisfies.
    """
    _require_keys(entry, _WITNESS_KEYS, "primal witness")
    if entry["objective_id"] != objective_id:
        raise CheckViolation(
            f"primal witness: objective {entry['objective_id']!r} is not "
            f"the record objective {objective_id!r}"
        )
    weights = entry["weights"]
    if not isinstance(weights, list) or not weights:
        raise CheckViolation("primal witness: carries no weights")
    objective = _route(objective_id)[0]
    total = balance = low_g = high_g = low_h = high_h = value = Fraction(0)
    previous = None
    for item in weights:
        _require_keys(item, _WEIGHT_KEYS, "witness weight")
        index = _exact_int(item["state_index"], "witness state_index")
        if not 0 <= index < len(states):
            raise CheckViolation(
                f"state_index: {index} is out of the state range"
            )
        if previous is not None and index <= previous:
            raise CheckViolation(
                f"state_index: {index} does not increase past {previous}"
            )
        previous = index
        weight = _canonical_fraction(
            item["weight"], f"weight of state {index}"
        )
        if weight < 0:
            raise CheckViolation(
                f"weight: {weight} of state {index} is negative"
            )
        state = states[index]
        total += weight
        balance += weight * state["excess_balance"]
        low_g += weight * state["g_lo"]
        high_g += weight * state["g_hi"]
        low_h += weight * state["h_lo"]
        high_h += weight * state["h_hi"]
        value += weight * objective(state)
    if total != VERTICES:
        raise CheckViolation(
            f"weight: witness weights sum to {total}, not {VERTICES}"
        )
    if balance != 0:
        raise CheckViolation(
            f"witness: aggregate excess balance is {balance}, not zero"
        )
    if low_g > 0:
        raise CheckViolation(
            f"witness: aggregate g_lo is {low_g}, which exceeds zero"
        )
    if high_g < 0:
        raise CheckViolation(
            f"witness: aggregate g_hi is {high_g}, which is below zero"
        )
    if low_h > 0:
        raise CheckViolation(
            f"witness: aggregate h_lo is {low_h}, which exceeds zero"
        )
    if high_h < 0:
        raise CheckViolation(
            f"witness: aggregate h_hi is {high_h}, which is below zero"
        )
    return value


def _stored_exact(text, expected, label):
    """A stored exact value exists exactly when its evidence does."""
    if expected is None:
        if text is not None:
            raise CheckViolation(f"{label}: {text!r} has no evidence")
        return
    if text is None:
        raise CheckViolation(
            f"{label}: is null but the evidence gives {expected}"
        )
    if _canonical_fraction(text, label) != expected:
        raise CheckViolation(
            f"{label}: {text!r} is not the verified {expected}"
        )


def _verify_route_record(states, record):
    """Verify one route record against the rebuilt states."""
    _require_keys(record, _SEARCH_RECORD_KEYS, "search record")
    objective_id = record["objective_id"]
    objective, limit, inclusive = _route(objective_id)
    bound = None
    if record["certificate"] is not None:
        bound = _verify_certificate(
            states, record["certificate"], objective_id
        )
    _stored_exact(record["exact_upper_bound"], bound, "exact_upper_bound")
    value = None
    if record["primal_witness"] is not None:
        value = _verify_witness(
            states, record["primal_witness"], objective_id
        )
    _stored_exact(record["exact_witness_value"], value, "exact_witness_value")
    status = record["route_status"]
    if status not in ROUTE_STATUSES:
        raise CheckViolation(
            f"route_status: {status!r} is not a route status"
        )
    accepted = bound is not None and _clears_edge(bound, limit, inclusive)
    rejected = value is not None and not _clears_edge(
        value, limit, inclusive
    )
    if accepted and rejected:
        raise CheckViolation(
            f"route_status: accepted bound {bound} contradicts the witness "
            f"value {value} of {objective_id}"
        )
    if accepted:
        derived = ACCEPTED
    elif rejected:
        derived = REJECTED
    else:
        derived = UNRESOLVED
    if status != derived:
        raise CheckViolation(
            f"route_status: {status!r} is not the derived {derived} of "
            f"{objective_id}"
        )
    return derived


def _check_searches(document, states):
    """Verify the search records against the rebuilt states.

    Returns the three route statuses, or None for the pre-search
    artifact (no searches yet).
    """
    searches = document["searches"]
    if not isinstance(searches, list):
        raise CheckViolation("searches: not a list")
    if not searches:
        return None                 # the pre-search variant
    expected = len(ROUTE_IDS) + 1
    if len(searches) != expected:
        raise CheckViolation(
            f"searches: {len(searches)} records, expected exactly {expected}"
        )
    statuses = []
    for index, objective_id in enumerate(ROUTE_IDS):
        record = searches[index]
        _require_keys(record, _SEARCH_RECORD_KEYS, "search record")
        if record["objective_id"] != objective_id:
            raise CheckViolation(
                f"route order: searches[{index}] objective "
                f"{record['objective_id']!r} is not the frozen route order "
                f"{objective_id!r}"
            )
        statuses.append(_verify_route_record(states, record))
    if searches[-1] != UNAVAILABLE_RECORD:
        raise CheckViolation(
            f"unavailable: searches[-1] must be the frozen "
            f"{UNAVAILABLE_ROUTE} record carrying no evidence, got "
            f"{searches[-1]!r}"
        )
    return statuses


def _terminal_disposition(statuses):
    """The verdict implied by the three route statuses alone."""
    if any(status == ACCEPTED for status in statuses):
        return DISPOSITION_ACCEPTED
    if all(status == REJECTED for status in statuses):
        return DISPOSITION_NO_CUT
    return DISPOSITION_UNRESOLVED


# ---------------------------------------------------------------- schema ---

_TOP_LEVEL_KEYS = (
    "campaign_id", "catalog_motif_histograms", "catalog_motif_windows",
    "disposition", "inputs", "n49_replay", "outer_motif_windows",
    "r35_edge_windows", "replay", "schema_version", "searches", "states",
    "trust_roots",
)
_INPUT_KEYS = ("bytes", "graph_count", "relative_path", "sha256")
_TRUST_KEYS = ("statement", "url")
_REPLAY_KEYS = ("complements", "published_graphs", "residual_zero")
_N49_KEYS = ("extremal_graphs", "forced_i4_mean", "observed_i4", "row_scaled")
_R35_WINDOW_KEYS = ("edge_max", "edge_min", "order")
_HISTOGRAM_KEYS = ("edges", "graph_count", "histograms", "order", "source")
_MOTIF_WINDOW_KEYS = ("edges", "order", "source", "windows")
_OUTER_WINDOW_KEYS = ("edges", "method", "order", "windows")
_HISTOGRAM_ENTRY_KEYS = ("count", "value")
_STATE_KEYS = (
    "a", "b", "d", "deficiency", "excess_balance", "g_hi", "g_lo", "h_hi",
    "h_lo", "x_interval_source", "x_motif_source", "y_interval_source",
    "y_motif_source",
)
_SEARCH_RECORD_KEYS = (
    "certificate", "exact_upper_bound", "exact_witness_value", "objective_id",
    "primal_witness", "route_status",
)
_CERTIFICATE_KEYS = (
    "alpha", "beta", "delta", "epsilon", "gamma", "objective_id", "zeta",
)
_WITNESS_KEYS = ("objective_id", "weights")
_WEIGHT_KEYS = ("state_index", "weight")


def _check_relative_path(record):
    """A canonical `r55/...` POSIX path; nothing absolute may appear."""
    label = record["relative_path"]
    if type(label) is not str or not label:
        raise CheckViolation("relative_path: input record without a path")
    pure = PurePosixPath(label)
    if pure.is_absolute() or label.startswith("/"):
        raise CheckViolation(f"relative_path: {label!r} is absolute")
    parts = pure.parts
    if not parts or parts[0] != "r55":
        raise CheckViolation(
            f"relative_path: {label!r} is not rooted at r55/"
        )
    if any(part in ("", ".", "..") for part in parts):
        raise CheckViolation(f"relative_path: {label!r} is not canonical")
    return label


def _check_window_pair(record, motif, label, cap):
    """One [lo, hi] motif window: exact ints, ordered, inside the cap."""
    pair = record["windows"][motif]
    if not isinstance(pair, list) or len(pair) != 2:
        raise CheckViolation(
            f"{label}: window {motif} is not a [lo, hi] pair"
        )
    low = _exact_int(pair[0], f"{label} {motif} lo")
    high = _exact_int(pair[1], f"{label} {motif} hi")
    if low < 0 or low > high:
        raise CheckViolation(
            f"{label}: window {motif} is empty or negative "
            f"([{low}, {high}])"
        )
    if high > cap:
        raise CheckViolation(
            f"{label}: window {motif} hi {high} exceeds the sound bound "
            f"{cap}"
        )
    return low, high


def _check_schema(document):
    """Exact key sets and per-section self-consistency, no disk access."""
    if not isinstance(document, dict):
        raise CheckViolation(f"schema: document is not an object")
    _require_keys(document, _TOP_LEVEL_KEYS, "schema: analysis document")
    if _exact_int(document["schema_version"], "schema_version") != (
        SCHEMA_VERSION
    ):
        raise CheckViolation(
            f"schema_version: {document['schema_version']!r}"
        )
    if document["campaign_id"] != CAMPAIGN_ID:
        raise CheckViolation(f"campaign_id: {document['campaign_id']!r}")
    if type(document["disposition"]) is not str:
        raise CheckViolation(
            f"disposition: {document['disposition']!r} is not a string"
        )
    for text in _iter_strings(document):
        if text.startswith("/"):
            raise CheckViolation(
                f"schema: absolute path {text!r} in the analysis"
            )

    inputs = document["inputs"]
    if not isinstance(inputs, list) or not inputs:
        raise CheckViolation("inputs: no input records")
    paths = []
    for record in inputs:
        _require_keys(record, _INPUT_KEYS, "input record")
        if not _is_sha256(record["sha256"]):
            raise CheckViolation(
                f"sha256: {record.get('relative_path')!r} is malformed"
            )
        if _exact_int(record["bytes"], "bytes") <= 0:
            raise CheckViolation("bytes: empty input")
        count = record["graph_count"]
        if count is not None and _exact_int(count, "graph_count") <= 0:
            raise CheckViolation("graph_count: non-positive")
        paths.append(_check_relative_path(record))
    if paths != sorted(set(paths)):
        raise CheckViolation(
            "inputs: records are not sorted and unique by relative_path"
        )

    roots = document["trust_roots"]
    if not isinstance(roots, list) or len(roots) != TRUST_ROOT_COUNT:
        raise CheckViolation(
            f"trust_roots: expected exactly {TRUST_ROOT_COUNT} statements"
        )
    for record in roots:
        _require_keys(record, _TRUST_KEYS, "trust root")
        if type(record["statement"]) is not str:
            raise CheckViolation(
                f"trust_roots: statement {record['statement']!r} is not "
                "a string"
            )
        if not record["statement"]:
            raise CheckViolation("trust_roots: root without a statement")
        if type(record["url"]) is not str:
            raise CheckViolation(
                f"trust_roots: url {record['url']!r} is not a string"
            )
        if not record["url"].startswith("https://"):
            raise CheckViolation(
                f"trust_roots: url {record['url']!r} is not https"
            )

    replay = document["replay"]
    _require_keys(replay, _REPLAY_KEYS, "replay")
    for key in _REPLAY_KEYS:
        _exact_int(replay[key], f"replay {key}")

    n49 = document["n49_replay"]
    _require_keys(n49, _N49_KEYS, "n49_replay")
    extremal_count = _exact_int(n49["extremal_graphs"], "extremal_graphs")
    if extremal_count <= 0:
        raise CheckViolation("extremal_graphs: non-positive")
    _exact_int(n49["row_scaled"], "row_scaled")
    _exact_int(n49["forced_i4_mean"], "forced_i4_mean")
    observed = n49["observed_i4"]
    if not isinstance(observed, list) or not observed:
        raise CheckViolation("observed_i4: no values")
    for value in observed:
        _exact_int(value, "observed_i4 value")

    r35_windows = document["r35_edge_windows"]
    if not isinstance(r35_windows, list):
        raise CheckViolation("r35_edge_windows: not a list")
    orders = []
    for record in r35_windows:
        _require_keys(record, _R35_WINDOW_KEYS, "R(3,5) window")
        order = _exact_int(record["order"], "R(3,5) window order")
        low = _exact_int(record["edge_min"], "R(3,5) window edge_min")
        high = _exact_int(record["edge_max"], "R(3,5) window edge_max")
        if low > high:
            raise CheckViolation(
                f"r35_edge_windows: order {order} window is empty"
            )
        orders.append(order)
    if orders != list(range(R35_MAX_ORDER + 1)):
        raise CheckViolation(
            f"r35_edge_windows: orders {orders} do not cover exactly "
            f"0..{R35_MAX_ORDER}"
        )

    histograms = document["catalog_motif_histograms"]
    if not isinstance(histograms, list) or not histograms:
        raise CheckViolation("catalog_motif_histograms: no classes")
    if len(histograms) <= R45_EXTREME_RECORDS:
        raise CheckViolation(
            f"catalog_motif_histograms: {len(histograms)} classes, expected "
            f"more than the {R45_EXTREME_RECORDS} extreme classes"
        )
    observed_windows = {}
    classes = []
    sources = set()
    for record in histograms:
        _require_keys(record, _HISTOGRAM_KEYS, "catalog histogram")
        order = _exact_int(record["order"], "catalog histogram order")
        edges = _exact_int(record["edges"], "catalog histogram edges")
        key = (order, edges)
        classes.append(key)
        source = record["source"]
        if not isinstance(source, str) or not source.endswith(".g6"):
            raise CheckViolation(
                f"source: class {key} has no graph6 source"
            )
        sources.add(source)
        section = record["histograms"]
        if (not isinstance(section, dict)
                or tuple(sorted(section)) != tuple(sorted(MOTIFS))):
            raise CheckViolation(
                f"catalog_motif_histograms: class {key} does not cover "
                f"exactly the five motifs"
            )
        graph_count = _exact_int(
            record["graph_count"], "catalog histogram graph_count"
        )
        if graph_count <= 0:
            raise CheckViolation(
                f"graph_count: class {key} is non-positive"
            )
        span = {}
        for motif in MOTIFS:
            entries = section[motif]
            if not isinstance(entries, list) or not entries:
                raise CheckViolation(
                    f"catalog_motif_histograms: class {key} has an empty "
                    f"{motif} histogram"
                )
            values = []
            total = 0
            for entry in entries:
                _require_keys(entry, _HISTOGRAM_ENTRY_KEYS, "motif count")
                frequency = _exact_int(
                    entry["count"], f"{motif} histogram count"
                )
                if frequency <= 0:
                    raise CheckViolation(
                        f"catalog_motif_histograms: class {key} has a "
                        "non-positive frequency"
                    )
                values.append(_exact_int(entry["value"], f"{motif} value"))
                total += frequency
            if values != sorted(set(values)):
                raise CheckViolation(
                    f"catalog_motif_histograms: class {key} has unsorted "
                    f"or repeated {motif} values"
                )
            if total != graph_count:
                raise CheckViolation(
                    f"graph_count: class {key} {motif} histogram sums to "
                    f"{total}, not its {graph_count} graphs"
                )
            span[motif] = (values[0], values[-1])
        observed_windows[key] = (span, source)
    if classes != sorted(set(classes)):
        raise CheckViolation(
            "catalog_motif_histograms: classes are not sorted and unique"
        )

    exact_windows = document["catalog_motif_windows"]
    if not isinstance(exact_windows, list):
        raise CheckViolation("catalog_motif_windows: not a list")
    if len(exact_windows) != len(classes):
        raise CheckViolation(
            f"catalog_motif_windows: {len(exact_windows)} windows for "
            f"{len(classes)} classes"
        )
    seen = []
    for record in exact_windows:
        _require_keys(record, _MOTIF_WINDOW_KEYS, "catalog window")
        order = _exact_int(record["order"], "catalog window order")
        edges = _exact_int(record["edges"], "catalog window edges")
        key = (order, edges)
        if key not in observed_windows:
            raise CheckViolation(
                f"catalog_motif_windows: window {key} has no histogram"
            )
        span, histogram_source = observed_windows[key]
        windows = record["windows"]
        if (not isinstance(windows, dict)
                or tuple(sorted(windows)) != tuple(sorted(MOTIFS))):
            raise CheckViolation(
                f"catalog_motif_windows: window {key} does not cover "
                f"exactly the five motifs"
            )
        for motif in MOTIFS:
            pair = windows[motif]
            if not isinstance(pair, list) or len(pair) != 2:
                raise CheckViolation(
                    f"catalog_motif_windows: window {key} motif {motif} "
                    "is not a [lo, hi] pair"
                )
            low = _exact_int(pair[0], f"catalog window {motif} lo")
            high = _exact_int(pair[1], f"catalog window {motif} hi")
            if (low, high) != span[motif]:
                raise CheckViolation(
                    f"catalog_motif_windows: window {key} motif {motif} "
                    f"[{low}, {high}] disagrees with its own histogram "
                    f"{span[motif]}"
                )
        if record["source"] != histogram_source:
            raise CheckViolation(
                f"source: catalog window {key} records "
                f"{record['source']!r}, its histogram records "
                f"{histogram_source!r}"
            )
        seen.append(key)
    if seen != sorted(set(seen)):
        raise CheckViolation(
            "catalog_motif_windows: classes are not sorted and unique"
        )

    outer = document["outer_motif_windows"]
    if not isinstance(outer, list) or not outer:
        raise CheckViolation("outer_motif_windows: no classes")
    outer_classes = []
    for record in outer:
        _require_keys(record, _OUTER_WINDOW_KEYS, "outer window")
        order = _exact_int(record["order"], "outer window order")
        edges = _exact_int(record["edges"], "outer window edges")
        key = (order, edges)
        if record["method"] not in OUTER_METHODS:
            raise CheckViolation(
                f"outer_motif_windows: class {key} carries the unknown "
                f"method {record['method']!r}"
            )
        windows = record["windows"]
        if (not isinstance(windows, dict)
                or tuple(sorted(windows)) != tuple(sorted(MOTIFS))):
            raise CheckViolation(
                f"outer_motif_windows: window {key} does not cover exactly "
                f"the five motifs"
            )
        for motif in MOTIFS:
            _check_window_pair(
                record, motif, "outer_motif_windows",
                _trivial_cap(order, motif),
            )
        outer_classes.append(key)
    if outer_classes != sorted(set(outer_classes)):
        raise CheckViolation(
            "outer_motif_windows: classes are not sorted and unique"
        )

    states = document["states"]
    if not isinstance(states, list) or not states:
        raise CheckViolation("states: the analysis records no states")
    cells = []
    for record in states:
        _require_keys(record, _STATE_KEYS, "state")
        for field in ("a", "b", "d", "deficiency", "excess_balance",
                      "g_hi", "g_lo", "h_hi", "h_lo"):
            _exact_int(record[field], f"state {field}")
        if record["d"] not in STATE_DEGREES:
            raise CheckViolation(
                f"d: state degree {record['d']} outside "
                f"{STATE_DEGREES[0]}..{STATE_DEGREES[-1]}"
            )
        if record["a"] < 0 or record["b"] < 0:
            raise CheckViolation(
                f"state: negative shortfall at d={record['d']}"
            )
        if record["deficiency"] != record["a"] + record["b"]:
            raise CheckViolation(
                "deficiency: disagrees with a + b"
            )
        if record["g_lo"] > record["g_hi"]:
            raise CheckViolation("g_lo: state has an empty g interval")
        if record["h_lo"] > record["h_hi"]:
            raise CheckViolation("h_lo: state has an empty h interval")
        for field in ("x_interval_source", "y_interval_source",
                      "x_motif_source", "y_motif_source"):
            if type(record[field]) is not str or not record[field]:
                raise CheckViolation(
                    f"{field}: empty provenance string"
                )
        cells.append((record["d"], record["a"], record["b"]))
    if cells != sorted(set(cells)):
        raise CheckViolation("states: not sorted and unique")

    searches = document["searches"]
    if not isinstance(searches, list):
        raise CheckViolation("searches: not a list")
    if searches and len(searches) != len(ROUTE_IDS) + 1:
        raise CheckViolation(
            f"searches: expected exactly {len(ROUTE_IDS) + 1} records"
        )


# ---------------------------------------------------------------- hashes ---

def _check_hashes(document, root):
    """Every recorded input still hashes to its recorded SHA-256."""
    for record in document["inputs"]:
        label = record["relative_path"]
        path = root / label
        if not path.is_file():
            raise CheckViolation(
                f"inputs: {label} is missing under the root {root}"
            )
        digest = _sha256_of(path)
        if digest != record["sha256"]:
            raise CheckViolation(
                f"sha256: {label} records {record['sha256']}, the local "
                f"bytes hash to {digest}"
            )


def _check_sizes(document, root):
    """Every recorded byte count equals the file size on disk.

    Redundant with the digest in isolation, but it pins the recorded
    `bytes` field itself; it runs after the validation selection so a
    structurally broken report is named before any byte count is.
    """
    for record in document["inputs"]:
        label = record["relative_path"]
        try:
            size = (root / label).stat().st_size
        except OSError as exc:
            raise CheckViolation(
                f"bytes: {label}: unstatable input ({exc})"
            ) from None
        if size != record["bytes"]:
            raise CheckViolation(
                f"bytes: {label} records {record['bytes']}, the local file "
                f"holds {size}"
            )


# ---------------------------------------------------------- catalog stage ---

def _derived_inputs(main_records, extreme_records):
    """The required input set: rel path -> required graph_count or None."""
    required = {}
    for order, name in enumerate(R35_FILES, start=1):
        required[f"r55/{DATA_DIRNAME}/{name}"] = (
            main_records[name]["graphs"]
        )
    for name in sorted(extreme_records):
        required[
            f"r55/{DATA_DIRNAME}/{R45_EXTREME_DIRNAME}/{name}"
        ] = extreme_records[name]["graphs"]
    required[f"r55/{DATA_DIRNAME}/{R45_CENSUS_FILE}"] = (
        main_records[R45_CENSUS_FILE]["graphs"]
    )
    required[f"r55/{DATA_DIRNAME}/{R55_REPLAY_FILE}"] = (
        main_records[R55_REPLAY_FILE]["graphs"]
    )
    for name in (VALIDATION_FILE, VALIDATION_EXTREME_FILE, TABLE_FILE):
        required[f"r55/{DATA_DIRNAME}/{name}"] = None
    for name in RECORDED_SOURCES:
        required[f"r55/src/{name}"] = None
    return required


def _check_inputs_vs_derived(document, required):
    """The recorded inputs are exactly the derived required set.

    The frozen m=3 machinery that the producer imports, the shared
    identity module, and the read-only m=3 artifact may appear as
    additional hash-pinned provenance; anything else is unexpected.
    """
    recorded = {
        record["relative_path"]: record["graph_count"]
        for record in document["inputs"]
    }
    missing = sorted(set(required) - set(recorded))
    if missing:
        raise CheckViolation(f"inputs: no record for {missing}")
    extra = sorted(set(recorded) - set(required) - set(OPTIONAL_INPUTS))
    if extra:
        raise CheckViolation(f"inputs: unexpected records {extra}")
    for label in sorted(required):
        if recorded[label] != required[label]:
            raise CheckViolation(
                f"graph_count: {label} records {recorded[label]}, the "
                f"validated files give {required[label]}"
            )
    for label in OPTIONAL_INPUTS:
        if label in recorded and recorded[label] is not None:
            raise CheckViolation(
                f"graph_count: {label} records {recorded[label]}, "
                "provenance-only inputs carry no graph count"
            )


def _check_catalogs(document, root):
    """Re-select, re-count, and re-sweep every catalog the artifact claims.

    Returns the derived (windows, owners, catalog q windows, motif
    windows, streamed graph total) for the state rebuild.
    """
    data = root / "r55" / DATA_DIRNAME
    main_records = select_main_records(data / VALIDATION_FILE)
    windows = select_r35_windows(data / VALIDATION_FILE)
    extreme_records = select_extreme_records(data / VALIDATION_EXTREME_FILE)

    _check_inputs_vs_derived(
        document, _derived_inputs(main_records, extreme_records)
    )

    _check_sizes(document, root)

    for order, name in enumerate(R35_FILES, start=1):
        graphs = count_graph6_lines(data / name)
        if graphs != main_records[name]["graphs"]:
            raise CheckViolation(
                f"count: {name} holds {graphs} non-empty lines, the "
                f"validation report records "
                f"{main_records[name]['graphs']}"
            )

    replay_lines = count_graph6_lines(data / R55_REPLAY_FILE)
    if replay_lines != main_records[R55_REPLAY_FILE]["graphs"]:
        raise CheckViolation(
            f"count: {R55_REPLAY_FILE} holds {replay_lines} non-empty "
            f"lines, the validation report records "
            f"{main_records[R55_REPLAY_FILE]['graphs']}"
        )

    recorded_shas = {
        record["relative_path"]: record["sha256"]
        for record in document["inputs"]
    }
    for name in (*R35_FILES, R45_CENSUS_FILE, R55_REPLAY_FILE):
        label = f"r55/{DATA_DIRNAME}/{name}"
        validated = main_records[name]["sha256"]
        if recorded_shas[label] != validated:
            raise CheckViolation(
                f"sha256: {label}: the artifact records "
                f"{recorded_shas[label]}, the validation report records "
                f"{validated}"
            )

    records = []
    for name in sorted(extreme_records):
        info = extreme_records[name]
        records.append((
            data / R45_EXTREME_DIRNAME / name,
            info["order"], info["edges"], info["graphs"],
            f"r55/{DATA_DIRNAME}/{R45_EXTREME_DIRNAME}/{name}",
        ))
    census_rel = f"r55/{DATA_DIRNAME}/{R45_CENSUS_FILE}"
    records.append((
        data / R45_CENSUS_FILE,
        main_records[R45_CENSUS_FILE]["order"], None,
        main_records[R45_CENSUS_FILE]["graphs"], census_rel,
    ))
    histograms, owners, _counts = sweep_catalogs(records)

    for key, owner in sorted(owners.items()):
        order, edges = key
        if owner == census_rel:
            if order != CENSUS_ORDER or not (
                CENSUS_EDGE_BAND[0] <= edges <= CENSUS_EDGE_BAND[1]
            ):
                raise CheckViolation(
                    f"owner: the census owns (order={order}, edges={edges}), "
                    f"outside the order-{CENSUS_ORDER} band "
                    f"{CENSUS_EDGE_BAND}"
                )
        elif order == CENSUS_ORDER:
            raise CheckViolation(
                f"owner: order-{CENSUS_ORDER} class {key} is owned by "
                f"{owner}, not the census"
            )

    doc_histograms = document["catalog_motif_histograms"]
    doc_classes = [(r["order"], r["edges"]) for r in doc_histograms]
    observed_classes = sorted(histograms)
    if doc_classes != observed_classes:
        raise CheckViolation(
            f"catalog_motif_histograms: recorded classes {doc_classes} are "
            f"not the observed {observed_classes}"
        )
    for record in doc_histograms:
        key = (record["order"], record["edges"])
        frequency = histograms[key]
        expected = {
            motif: [
                {"value": value, "count": frequency[motif][value]}
                for value in sorted(frequency[motif])
            ]
            for motif in MOTIFS
        }
        if record["histograms"] != expected:
            raise CheckViolation(
                f"catalog_motif_histograms: class (order={key[0]}, "
                f"edges={key[1]}) does not match the reswept motif "
                f"frequencies; recorded {record['histograms']}, observed "
                f"{expected}"
            )
        graphs_in_class = sum(frequency["q"].values())
        if record["graph_count"] != graphs_in_class:
            raise CheckViolation(
                f"graph_count: class (order={key[0]}, edges={key[1]}) "
                f"records {record['graph_count']}, the resweep gives "
                f"{graphs_in_class}"
            )
        if record["source"] != owners[key]:
            raise CheckViolation(
                f"source: class (order={key[0]}, edges={key[1]}) records "
                f"{record['source']!r}, the observed owner is "
                f"{owners[key]!r}"
            )

    for record in document["catalog_motif_windows"]:
        key = (record["order"], record["edges"])
        frequency = histograms[key]
        for motif in MOTIFS:
            pair = record["windows"][motif]
            if pair != [min(frequency[motif]), max(frequency[motif])]:
                raise CheckViolation(
                    f"catalog_motif_windows: class (order={key[0]}, "
                    f"edges={key[1]}) records {motif} window {pair}, "
                    f"observed [{min(frequency[motif])}, "
                    f"{max(frequency[motif])}]"
                )
        if record["source"] != owners[key]:
            raise CheckViolation(
                f"source: class (order={key[0]}, edges={key[1]}) window "
                f"source is not the observed owner {owners[key]!r}"
            )

    for record in document["outer_motif_windows"]:
        key = (record["order"], record["edges"])
        frequency = histograms.get(key)
        if frequency is None:
            continue        # a missing stratum: no data to audit against
        for motif in MOTIFS:
            low, high = record["windows"][motif]
            outside = [
                value for value in frequency[motif]
                if not low <= value <= high
            ]
            if outside:
                raise CheckViolation(
                    f"outer_motif_windows: class (order={key[0]}, "
                    f"edges={key[1]}) records {motif} window "
                    f"[{low}, {high}] excluding the observed value(s) "
                    f"{outside}"
                )

    for key in observed_classes:
        try:
            low, high = q_outer_interval(key[0], key[1], windows)
        except InfeasibleEdgeClass as exc:
            raise CheckViolation(
                f"envelope: class (order={key[0]}, edges={key[1]}) has no "
                f"surviving degree histogram ({exc})"
            ) from None
        outside = [
            q for q in histograms[key]["q"] if not low <= q <= high
        ]
        if outside:
            raise CheckViolation(
                f"envelope: class (order={key[0]}, edges={key[1]}) observes "
                f"q values {outside} outside its outer interval "
                f"[{low}, {high}]"
            )

    catalog_windows = {
        key: (min(histograms[key]["q"]), max(histograms[key]["q"]))
        for key in observed_classes
    }
    motif_windows = {
        key: {
            motif: (min(histograms[key][motif]), max(histograms[key][motif]))
            for motif in MOTIFS
        }
        for key in observed_classes
    }
    total = sum(
        sum(frequency["q"].values()) for frequency in histograms.values()
    )
    return windows, owners, catalog_windows, motif_windows, total


# ------------------------------------------------------------ state stage ---

def _check_states(document, root, windows, owners, catalog_windows,
                  motif_windows):
    """Rebuild the whole cone and require every recorded field to match."""
    edge_bounds = load_edge_bounds(
        root / "r55" / DATA_DIRNAME / TABLE_FILE
    )

    derived_windows = [
        {
            "order": order,
            "edge_min": windows[order][0],
            "edge_max": windows[order][1],
        }
        for order in sorted(windows)
    ]
    if document["r35_edge_windows"] != derived_windows:
        raise CheckViolation(
            "r35_edge_windows: the recorded windows are not the validated "
            f"ones ({derived_windows})"
        )

    outer_windows = {
        (record["order"], record["edges"]): {
            motif: (record["windows"][motif][0], record["windows"][motif][1])
            for motif in MOTIFS
        }
        for record in document["outer_motif_windows"]
    }
    rebuilt = rebuild_states(
        windows, edge_bounds, catalog_windows, owners, motif_windows,
        outer_windows
    )
    states = document["states"]
    if len(states) != len(rebuilt):
        raise CheckViolation(
            f"states: {len(states)} recorded, the rebuild gives "
            f"{len(rebuilt)}"
        )
    for index, (record, mine) in enumerate(zip(states, rebuilt)):
        label = (
            f"states[{index}] (d={record['d']}, a={record['a']}, "
            f"b={record['b']})"
        )
        for field in ("d", "a", "b", "deficiency", "excess_balance",
                      "g_lo", "g_hi", "h_lo", "h_hi"):
            if record[field] != mine[field]:
                raise CheckViolation(
                    f"{field}: {label} records {record[field]!r}, the "
                    f"rebuild gives {mine[field]!r}"
                )
        for field in ("x_interval_source", "y_interval_source"):
            if record[field] != mine[field]:
                raise CheckViolation(
                    f"{field}: {label} records {record[field]!r}, the "
                    f"rebuild gives {mine[field]!r}"
                )
        for field, acceptable in (
            ("x_motif_source", mine["x_motif_ok"]),
            ("y_motif_source", mine["y_motif_ok"]),
        ):
            if record[field] not in acceptable:
                raise CheckViolation(
                    f"{field}: {label} records {record[field]!r}, expected "
                    f"one of {sorted(acceptable)}"
                )
    return rebuilt


# ---------------------------------------------------------------- verdict ---

def verify_analysis(document, root=None, echo=None):
    """Re-verify one analysis document end to end; return the disposition.

    Every stage fails fast with a `CheckViolation` naming the first
    broken field.  `root` substitutes the math repository root that all
    recorded `r55/...` paths resolve against; it defaults to the
    checker's own repository root.  `echo`, when given, receives one
    deterministic progress line per completed stage.
    """
    progress = echo if callable(echo) else (lambda line: None)
    if root is None:
        root = _DEFAULT_ROOT
    root = Path(root).resolve()

    _check_schema(document)
    _check_hashes(document, root)
    progress(f"M4 HASHES: {len(document['inputs'])} inputs byte-identical")

    windows, owners, catalog_windows, motif_windows, total = _check_catalogs(
        document, root
    )
    progress(
        f"M4 CATALOGS: {total} graphs reswept into {len(owners)} classes"
    )

    _check_n49(document, root)
    n49 = document["n49_replay"]
    progress(
        f"M4 N49 REPLAY: extremal={n49['extremal_graphs']} "
        f"row={n49['row_scaled']} forced_i4_mean={n49['forced_i4_mean']} "
        f"observed={n49['observed_i4']}"
    )

    states = _check_states(
        document, root, windows, owners, catalog_windows, motif_windows
    )
    progress(f"M4 STATES: {len(states)} states rebuilt exactly")

    statuses = _check_searches(document, states)
    if statuses is None:
        progress("M4 SEARCHES: pre-search artifact, no routes recorded")
    else:
        progress(
            f"M4 SEARCHES: {len(ROUTE_IDS)} frozen routes verified against "
            f"{len(states)} rebuilt states"
        )

    lines, zeros = _check_replay(document, root)
    progress(f"M4 REPLAY: {zeros} of {2 * lines} residuals zero")

    if statuses is None:
        disposition = DISPOSITION_ANALYZED
    else:
        disposition = _terminal_disposition(statuses)
    if document["disposition"] != disposition:
        raise CheckViolation(
            f"disposition: {document['disposition']!r} is not the terminal "
            f"{disposition!r}"
        )
    _check_hashes(document, root)
    return disposition


def main(argv=None):
    """Verify one analysis artifact; print the verdict, return the code."""
    parser = argparse.ArgumentParser(
        description=(
            "Independently re-verify an m=4 higher-identity analysis "
            "artifact against its recorded inputs."
        )
    )
    parser.add_argument(
        "--root",
        default=None,
        help=(
            "math repository root every recorded r55/... path resolves "
            "against (default: the checker's own repository root)"
        ),
    )
    parser.add_argument(
        "analysis",
        help="path to higher_identity_m4.json",
    )
    args = parser.parse_args(argv)
    root = Path(args.root) if args.root is not None else _DEFAULT_ROOT
    try:
        target = Path(args.analysis)
        try:
            document = json.loads(target.read_text(encoding="ascii"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise CheckViolation(
                f"analysis: {target}: unreadable JSON ({exc})"
            ) from None
        disposition = verify_analysis(document, root, echo=_echo)
    except CheckViolation as exc:
        print(f"M4 EVIDENCE FAILED: {exc}", file=sys.stderr)
        return 1
    except (OSError, ValueError) as exc:
        print(f"M4 EVIDENCE FAILED: analysis: {exc}", file=sys.stderr)
        return 1
    print(f"M4 EVIDENCE VERIFIED: {disposition}")
    return 0


def _echo(line):
    print(line, flush=True)


if __name__ == "__main__":
    sys.exit(main())
