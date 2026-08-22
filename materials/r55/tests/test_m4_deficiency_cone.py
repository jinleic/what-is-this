"""Failing envelope and state-cone tests for the m=4 campaign (Task 2 RED).

Task 2 Steps 1-2 of notes/higher_order_identity_m4_plan_2026-08-17.md:
every test below fails until math/r55/src/m4_deficiency_cone.py lands its
interfaces.  The twelve tests, in mandate order:

 1. Motif4VectorTests.test_motif4_vector_matches_naive_reference_on_catalogs
 2. AffineIntervalTests.test_affine_interval_sign_handling
 3. RelationSystemTests.test_relation_rows_match_plan_constants_and_hold_on_random_graphs
 4. WindowResolutionTests.test_catalog_windows_override_outer_windows
 5. StateConeTests.test_state_intervals_contain_exact_catalog_stratum_pairs
 6. WindowResolutionTests.test_outer_motif_windows_audit_every_observed_catalog_value
 7. StateConeTests.test_build_states_reproduces_frozen_m3_domain_and_records
 8. KernelReplayTests.test_kernel_replays_656_residuals_and_n49_ground_truth
 9. RelaxedEnvelopeTests.test_relaxed_rhs_intervals_equal_naive_histogram_aggregates
10. RelaxedEnvelopeTests.test_verified_basic_solution_count_is_independently_exact
11. RelaxedEnvelopeTests.test_enumeration_failure_falls_back_to_sound_trivial_window
12. StateConeTests.test_state_record_serializes_the_v2_schema

Pinned data policy for the self-loading cone interfaces
(`m4_state_interval(d, a, b)` / `build_states()`): the g-side reuses the
frozen m=3 machinery and must reproduce the 3,215 records of
data/higher_identity_m3.json exactly (7-field equality, test 7); the
h-side motif windows resolve in artifact order -- exact catalog windows
from data/higher_identity_m4.json when that artifact exists, otherwise
the sound trivial [0, C(order, 4)] windows.  Focused tests therefore run
in seconds without the ~0.8 h catalog stream and stay sound once the
artifact lands; tests assert h-soundness (h_lo <= exact h <= h_hi),
never artifact-dependent endpoints.

Whole-suite runtime target: <= 120 s.  Every clamp is noted at its use:

 * first 250 consecutive graphs of r4517.77.g6 plus the only two
   132-edge order-24 census graphs (single full scan of the 16 MB
   r45_24.g6, cached and shared by tests 1, 5, and 8);
 * 60 seeded random graphs on n <= 7 vertices;
 * <= 200 streamed graphs per class for the catalog-override windows;
 * <= 20 catalog graphs per stratum side, three fixture stratum pairs;
 * <= 300 streamed graphs per class for the outer-window audit;
 * one order-6 sweep for the relaxed right-hand sides;
 * one pinned class (5, 5) for the in-test C(17, 9) = 24,310-basis
   exact enumeration (the suite's largest single cost, ~5-15 s);
 * the frozen m=3 artifact is parsed once and cached module-side.
"""

import dataclasses
import itertools
import json
import math
import pathlib
import random
import sys
import unittest
from fractions import Fraction
from unittest import mock

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from check_ramsey import parse_graph6_line
from m3_deficiency_cone import (  # frozen, verified m=3 infrastructure
    degree_bounds, degree_histograms, load_edge_bounds, r35_edge_windows,
)
from m4_subgraph_identities import (  # the landed Task 1 kernel
    independent_quad_count, m4_residual_specialized, motif4,
)
from subgraph_identities import complement

import m4_deficiency_cone as m4
from m4_deficiency_cone import (  # expected to fail before Task 2 Step 3
    MOTIF_ORDER, State, affine_interval, build_states,
    catalog_motif_windows, motif_interval, motif4_vector,
    m4_state_interval, outer_motif_window, relaxed_rhs_intervals,
    relation_rows, state_record, verified_basic_solution_count,
)

DATA = ROOT / "data"
EXTREME = DATA / "r45extreme"
M3_ARTIFACT_PATH = DATA / "higher_identity_m3.json"
STATE_DEGREES = (20, 21, 22, 23, 24)
MOTIF4_KEYS = ("q", "t", "diamond", "k3k1", "i4")
STATE_FIELDS = ("d", "a", "b", "deficiency", "excess_balance",
                "g_lo", "g_hi", "h_lo", "h_hi")
# The plan's v2 state-record schema block enumerates exactly these keys
# (nine dataclass fields plus four provenance strings; the mandate's
# "12 keys" label is a miscount of the plan's own list).
STATE_RECORD_KEYS = (
    "d", "a", "b", "deficiency", "excess_balance", "g_lo", "g_hi",
    "h_lo", "h_hi", "x_interval_source", "y_interval_source",
    "x_motif_source", "y_motif_source",
)

# The eleven induced 4-vertex motifs in the plan's canonical order.
CANONICAL_MOTIFS = ("E4", "K2+2K1", "2K2", "P3+K1", "P4", "K1_3",
                    "K3+K1", "C4", "paw", "diamond", "K4")
# Canonical induced edge count, sum of squared motif degrees, and
# triangle multiplicity tau (K3+K1: 1, paw: 1, diamond: 2, K4: 4).
MOTIF_EDGES = {"E4": 0, "K2+2K1": 1, "2K2": 2, "P3+K1": 2, "P4": 3,
               "K1_3": 3, "K3+K1": 3, "C4": 4, "paw": 4, "diamond": 5,
               "K4": 6}
MOTIF_DEGSQ = {"E4": 0, "K2+2K1": 2, "2K2": 4, "P3+K1": 6, "P4": 10,
               "K1_3": 12, "K3+K1": 12, "C4": 16, "paw": 18,
               "diamond": 26, "K4": 36}
MOTIF_TAU = {"K3+K1": 1, "paw": 1, "diamond": 2, "K4": 4}

# (induced edges, descending degree sequence) -> motif name; these eleven
# pairs classify every 4-vertex graph.
_QUARTET_NAMES = {
    (0, (0, 0, 0, 0)): "E4",
    (1, (1, 1, 0, 0)): "K2+2K1",
    (2, (1, 1, 1, 1)): "2K2",
    (2, (2, 1, 1, 0)): "P3+K1",
    (3, (3, 1, 1, 1)): "K1_3",
    (3, (2, 2, 2, 0)): "K3+K1",
    (3, (2, 2, 1, 1)): "P4",
    (4, (2, 2, 2, 2)): "C4",
    (4, (3, 2, 2, 1)): "paw",
    (5, (3, 3, 2, 2)): "diamond",
    (6, (3, 3, 3, 3)): "K4",
}

# Present fixed-edge R(4,5) classes outside the cone strata, used for the
# catalog-window override (<= 200 graphs) and the outer-window audit
# (<= 300 graphs); streaming clamps noted at each use.
CATALOG_CLASSES = (
    (17, 41, "r4517.41.g6"),
    (18, 50, "r4518.50.g6"),
    (20, 68, "r4520.68.g6"),
)


def c2(n):
    return n * (n - 1) // 2


def c3(n):
    return n * (n - 1) * (n - 2) // 6


def c4(n):
    return n * (n - 1) * (n - 2) * (n - 3) // 24


def f_of(order, degree):
    """R4 right-side per-degree term f(d) = d*C(r-1-d,2) + 4*C(d,2)(r-1-d) + 9*C(d,3)."""
    return (degree * c2(order - 1 - degree)
            + 4 * c2(degree) * (order - 1 - degree)
            + 9 * c3(degree))


_SAMPLES = {}


def read_graphs(path, edges=None, limit=None):
    """Parsed graphs of one graph6 catalog, edge-filtered, test-cached.

    `limit` caps the retained (filtered) graphs; for the order-24 census
    the 132-edge filter matches only two graphs, so the 16 MB file is
    scanned once and the result is cached for tests 1, 5, and 8.
    """
    key = (str(path), edges, limit)
    if key in _SAMPLES:
        return _SAMPLES[key]
    graphs = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if not line.strip():
                continue
            _, adj = parse_graph6_line(line)
            if edges is not None:
                if sum(row.bit_count() for row in adj) // 2 != edges:
                    continue
            graphs.append(adj)
            if limit is not None and len(graphs) >= limit:
                break
    _SAMPLES[key] = graphs
    return graphs


_M3_STATES = None


def m3_artifact_states():
    """The 3,215 frozen m=3 state records, parsed once and cached."""
    global _M3_STATES
    if _M3_STATES is None:
        raw = json.loads(M3_ARTIFACT_PATH.read_text(encoding="ascii"))
        _M3_STATES = raw["states"]
    return _M3_STATES


# ------------------------------------------------- naive test references ---


def naive_quartet_counts(adj):
    """All eleven induced 4-vertex motif counts, by 4-subset enumeration."""
    n = len(adj)
    counts = dict.fromkeys(CANONICAL_MOTIFS, 0)
    for i, j, k, l in itertools.combinations(range(n), 4):
        aij = (adj[i] >> j) & 1
        aik = (adj[i] >> k) & 1
        ail = (adj[i] >> l) & 1
        ajk = (adj[j] >> k) & 1
        ajl = (adj[j] >> l) & 1
        akl = (adj[k] >> l) & 1
        di = aij + aik + ail
        dj = aij + ajk + ajl
        dk = aik + ajk + akl
        dl = ail + ajl + akl
        name = _QUARTET_NAMES[(aij + aik + ail + ajk + ajl + akl,
                               tuple(sorted((di, dj, dk, dl),
                                            reverse=True)))]
        counts[name] += 1
    return counts


def naive_triangles(adj):
    """Triangle count by triple enumeration."""
    n = len(adj)
    total = 0
    for a in range(n - 2):
        for b in range(a + 1, n - 1):
            if not (adj[a] >> b) & 1:
                continue
            for c in range(b + 1, n):
                if ((adj[a] >> c) & 1) and ((adj[b] >> c) & 1):
                    total += 1
    return total


def naive_motif4_vector(adj):
    """The five streaming aggregates, from the naive references only."""
    counts = naive_quartet_counts(adj)
    triangles = naive_triangles(adj)
    wedge = sum(c2(row.bit_count()) for row in adj)
    return {"q": wedge - triangles, "t": triangles,
            "diamond": counts["diamond"],
            "k3k1": counts["K3+K1"],
            "i4": counts["E4"]}


def random_graph(rng, order, density):
    adj = [0] * order
    for u, v in itertools.combinations(range(order), 2):
        if rng.random() < density:
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return adj


def relation_rhs(adj):
    """Exact R1-R5 right-hand sides of one graph, plan constants verbatim."""
    n = len(adj)
    degrees = [row.bit_count() for row in adj]
    edges = sum(degrees) // 2
    wedge = sum(c2(d) for d in degrees)
    return {
        "R1": c4(n),
        "R2": edges * c2(n - 2),
        "R3": edges * c2(n - 2) + 2 * wedge * (n - 4) + edges * (edges - 1),
        "R4": sum(f_of(n, d) for d in degrees),
        "R5": (n - 3) * naive_triangles(adj),
    }


def naive_relaxed_rhs(order, edges, windows):
    """Relaxed R2-R5 right-hand sides over the test's own histograms.

    Mirrors the imported m=3 q-machinery's feasible set exactly: every
    integer degree histogram of the class whose per-neighborhood triangle
    window is nonempty.  R2 is class-determined and stays a point; R3 is
    the W aggregate; R4 the f aggregate; R5 is t in [W - q_hi, W_hi - q_lo]
    with the class q-window endpoints, all min/max over surviving
    histograms.  Returns None when no histogram survives.
    """
    low, _high = degree_bounds(order)
    w_lo = w_hi = f_lo = f_hi = q_lo = q_hi = None
    for hist in degree_histograms(order, edges):
        wedge = tri3_lo = tri3_hi = fsum = 0
        for index, count in enumerate(hist):
            if not count:
                continue
            degree = low + index
            wedge += count * c2(degree)
            tri3_lo += count * windows[degree][0]
            tri3_hi += count * windows[degree][1]
            fsum += count * f_of(order, degree)
        tri_lo = -(-tri3_lo // 3)   # ceil; 3t >= tri3_lo, t integral
        tri_hi = tri3_hi // 3
        if tri_lo > tri_hi:
            continue
        w_lo = wedge if w_lo is None else min(w_lo, wedge)
        w_hi = wedge if w_hi is None else max(w_hi, wedge)
        f_lo = fsum if f_lo is None else min(f_lo, fsum)
        f_hi = fsum if f_hi is None else max(f_hi, fsum)
        q_lo = (wedge - tri_hi if q_lo is None
                else min(q_lo, wedge - tri_hi))
        q_hi = (wedge - tri_lo if q_hi is None
                else max(q_hi, wedge - tri_lo))
    if w_lo is None:
        return None
    exact_r2 = edges * c2(order - 2)
    return {"R2": (exact_r2, exact_r2),
            "R3": (w_lo, w_hi),
            "R4": (f_lo, f_hi),
            "R5": (w_lo - q_hi, w_hi - q_lo)}


def _solve_square_int(mat):
    """Exact solution of one square integer system, or None if singular.

    Fraction-free (Bareiss) forward elimination -- every intermediate
    entry stays an exact minor, so each division is exact -- followed by
    one Fraction back-substitution.
    """
    n = len(mat)
    a = [row[:] for row in mat]
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            replacement = None
            for i in range(k + 1, n):
                if a[i][k]:
                    replacement = i
                    break
            if replacement is None:
                return None
            a[k], a[replacement] = a[replacement], a[k]
        pivot = a[k][k]
        top = a[k]
        for i in range(k + 1, n):
            row = a[i]
            factor = row[k]
            for j in range(k + 1, n + 1):
                row[j] = (row[j] * pivot - top[j] * factor) // prev
            row[k] = 0
        prev = pivot
    x = [None] * n
    for i in range(n - 1, -1, -1):
        diag = a[i][i]
        if diag == 0:
            return None
        acc = Fraction(a[i][n])
        row = a[i]
        for j in range(i + 1, n):
            acc -= row[j] * x[j]
        x[i] = acc / diag
    return x


def enumerate_basic_solutions_exact(order, edges, rhs):
    """Independently enumerate every basic feasible solution of R1-R5.

    Vertices of {s >= 0 : R1, R2 exact, R3-R5 two-sided} come from
    activating 9 of the 17 inequality rows (11 nonnegativity plus 6 range
    endpoints): C(17, 9) = 24,310 candidate bases.  `rhs` holds the
    AGGREGATE intervals of naive_relaxed_rhs (W range, f range, t range);
    the two-sided ROW bounds below transform them into row space exactly
    as the relations dictate: R3 row = eq2 + 2W(r-4) + e(e-1), R4 row =
    the f aggregate itself, R5 row = (r-3)*t.  Each square integer
    system is solved fraction-free; feasibility rechecks s >= 0 and all
    three two-sided rows.  Returns (feasible_basis_count,
    {motif4 key: (exact min, exact max)} as Fractions), where t is the
    tau-row over (r-3), q is W(s) - t(s) with W from the R3 row, and
    diamond / k3k1 / i4 are raw coordinates.

    CLAMP: this in-test enumeration is the suite's largest single cost
    (~5-15 s); it runs once, for the one pinned small class (5, 5).
    """
    r1 = (1,) * 11
    r2 = tuple(MOTIF_EDGES[m] for m in CANONICAL_MOTIFS)
    r3 = tuple(MOTIF_EDGES[m] * MOTIF_EDGES[m] for m in CANONICAL_MOTIFS)
    r4 = tuple(MOTIF_DEGSQ[m] for m in CANONICAL_MOTIFS)
    r5 = tuple(MOTIF_TAU.get(m, 0) for m in CANONICAL_MOTIFS)
    eq1 = c4(order)
    eq2 = edges * c2(order - 2)
    w_lo, w_hi = rhs["R3"]
    t_lo, t_hi = rhs["R5"]
    ranges = (
        (r3, (eq2 + 2 * w_lo * (order - 4) + edges * (edges - 1),
              eq2 + 2 * w_hi * (order - 4) + edges * (edges - 1))),
        (r4, rhs["R4"]),
        (r5, ((order - 3) * t_lo, (order - 3) * t_hi)),
    )
    endpoints = tuple((k, side) for k in range(3) for side in (0, 1))
    pick_indices = (CANONICAL_MOTIFS.index("diamond"),
                    CANONICAL_MOTIFS.index("K3+K1"),
                    CANONICAL_MOTIFS.index("E4"))
    count = 0
    best = {key: None for key in MOTIF4_KEYS}

    def offer(key, value):
        current = best[key]
        if current is None:
            best[key] = (value, value)
        else:
            lo, hi = current
            best[key] = (lo if lo <= value else value,
                         hi if hi >= value else value)

    for zeroed_count in range(3, 10):
        active = 9 - zeroed_count
        for zeroed in itertools.combinations(range(11), zeroed_count):
            keep = [v for v in range(11) if v not in zeroed]
            eq1_row = [r1[v] for v in keep]
            eq1_row.append(eq1)
            eq2_row = [r2[v] for v in keep]
            eq2_row.append(eq2)
            restricted = [([row[v] for v in keep], bounds)
                          for row, bounds in ranges]
            pick_positions = tuple(
                None if i in zeroed else keep.index(i) for i in pick_indices)
            for chosen in itertools.combinations(endpoints, active):
                mat = [eq1_row[:], eq2_row[:]]
                for k, side in chosen:
                    row_k, bounds = restricted[k]
                    mat.append(row_k + [bounds[side]])
                sol = _solve_square_int(mat)
                if sol is None or any(value < 0 for value in sol):
                    continue
                feasible = True
                vals = []
                for row_k, bounds in restricted:
                    total = 0
                    for coef, x in zip(row_k, sol):
                        total += coef * x
                    if total < bounds[0] or total > bounds[1]:
                        feasible = False
                        break
                    vals.append(total)
                if not feasible:
                    continue
                count += 1
                picked = [Fraction(0) if pos is None else sol[pos]
                          for pos in pick_positions]
                offer("diamond", picked[0])
                offer("k3k1", picked[1])
                offer("i4", picked[2])
                t_val = vals[2] / (order - 3)
                offer("t", t_val)
                w_val = (vals[0] - eq2 - edges * (edges - 1)) \
                    / (2 * (order - 4))
                offer("q", w_val - t_val)
    return count, best


# ------------------------------------------------------------------ tests ---


class Motif4VectorTests(unittest.TestCase):
    def test_motif4_vector_matches_naive_reference_on_catalogs(self):
        # CLAMP: the first 250 consecutive graphs of the order-17 class
        # r4517.77 (naive 4-subset reference, ~2-3 s) plus the only two
        # 132-edge graphs of the order-24 census (single cached 16 MB
        # scan shared with tests 5 and 8).
        prefix = read_graphs(EXTREME / "r4517.77.g6", 77, 250)
        self.assertEqual(len(prefix), 250)
        for adj in prefix:
            vector = motif4_vector(adj)
            self.assertEqual(vector, naive_motif4_vector(adj))
            wedge = sum(c2(row.bit_count()) for row in adj)
            self.assertEqual(vector["q"], wedge - vector["t"])

        extremal = read_graphs(DATA / "r45_24.g6", 132, 20)
        self.assertEqual(len(extremal), 2)
        for adj in extremal:
            vector = motif4_vector(adj)
            self.assertEqual(vector, naive_motif4_vector(adj))
            self.assertEqual(vector["t"], 176)
            self.assertEqual(vector["diamond"], 792)
            self.assertEqual(vector["k3k1"], 528)
            self.assertIn(vector["i4"], (138, 144))
            wedge = sum(c2(row.bit_count()) for row in adj)
            self.assertEqual(vector["q"], wedge - vector["t"])
        self.assertEqual(
            sorted(motif4_vector(adj)["i4"] for adj in extremal),
            [138, 144],
        )


class AffineIntervalTests(unittest.TestCase):
    def test_affine_interval_sign_handling(self):
        self.assertEqual(affine_interval(-3, (10, 20)), (-60, -30))
        self.assertEqual(affine_interval(0, (10, 20)), (0, 0))
        self.assertEqual(affine_interval(2, (3, 5)), (6, 10))


class RelationSystemTests(unittest.TestCase):
    def test_relation_rows_match_plan_constants_and_hold_on_random_graphs(
            self):
        # The module's coefficient table is pinned motif-by-motif against
        # the plan's constants (edges, squared edges, degree-square sums,
        # triangle multiplicities), then R1-R5 are verified as exact
        # identities on 60 seeded random graphs, n <= 7, naive counts.
        self.assertEqual(tuple(MOTIF_ORDER), CANONICAL_MOTIFS)
        expected = {
            "R1": (1,) * 11,
            "R2": tuple(MOTIF_EDGES[m] for m in CANONICAL_MOTIFS),
            "R3": tuple(MOTIF_EDGES[m] ** 2 for m in CANONICAL_MOTIFS),
            "R4": tuple(MOTIF_DEGSQ[m] for m in CANONICAL_MOTIFS),
            "R5": tuple(MOTIF_TAU.get(m, 0) for m in CANONICAL_MOTIFS),
        }
        rows = relation_rows()
        self.assertEqual(rows, expected)

        # CLAMP: 60 seeded random graphs (seed 20260817), n <= 7.
        rng = random.Random(20260817)
        checked = 0
        for _ in range(60):
            order = rng.choice([5, 6, 7])
            density = rng.choice([0.3, 0.5, 0.7])
            adj = random_graph(rng, order, density)
            counts = naive_quartet_counts(adj)
            svec = [counts[m] for m in CANONICAL_MOTIFS]
            rhs = relation_rhs(adj)
            for name in ("R1", "R2", "R3", "R4", "R5"):
                self.assertEqual(
                    sum(c * s for c, s in zip(rows[name], svec)),
                    rhs[name], (name, adj),
                )
                checked += 1
        self.assertEqual(checked, 300)


class WindowResolutionTests(unittest.TestCase):
    def test_catalog_windows_override_outer_windows(self):
        # CLAMP: <= 200 streamed graphs per class, three present classes.
        for order, edges, name in CATALOG_CLASSES:
            graphs = read_graphs(EXTREME / name, edges, 200)
            self.assertTrue(graphs)
            vectors = [motif4_vector(adj) for adj in graphs]
            histograms = {key: {} for key in MOTIF4_KEYS}
            for vector in vectors:
                for key in MOTIF4_KEYS:
                    value = vector[key]
                    histograms[key][value] = \
                        histograms[key].get(value, 0) + 1
            records = {(order, edges): histograms}
            catalog = catalog_motif_windows(records)
            self.assertEqual(sorted(catalog), [(order, edges)])
            windows = catalog[(order, edges)]
            self.assertEqual(sorted(windows), sorted(MOTIF4_KEYS))
            stricter = 0
            for key in MOTIF4_KEYS:
                observed = (min(v[key] for v in vectors),
                            max(v[key] for v in vectors))
                self.assertEqual(windows[key], observed)
                # A present class resolves to its exact catalog window.
                self.assertEqual(
                    motif_interval(order, edges, key, catalog), observed)
                # A missing class resolves to the outer LP window.
                outer = motif_interval(order, edges, key, None)
                self.assertEqual(outer, outer_motif_window(order, edges, key))
                self.assertEqual(
                    motif_interval(order, edges, key, {}), outer)
                # ... which soundly contains the catalog window.
                self.assertLessEqual(outer[0], observed[0])
                self.assertLessEqual(observed[1], outer[1])
                if outer != observed:
                    stricter += 1
            # The relaxation is genuinely wider somewhere in every class.
            self.assertGreater(stricter, 0, (order, edges))

    def test_outer_motif_windows_audit_every_observed_catalog_value(self):
        # Audit direction on present classes: every observed motif value
        # of the streamed prefix lies inside its outer LP window.
        # CLAMP: <= 300 streamed graphs per class, the same three classes.
        for order, edges, name in CATALOG_CLASSES:
            graphs = read_graphs(EXTREME / name, edges, 300)
            self.assertTrue(graphs)
            vectors = [motif4_vector(adj) for adj in graphs]
            for key in MOTIF4_KEYS:
                lo, hi = outer_motif_window(order, edges, key)
                self.assertIs(type(lo), int)
                self.assertIs(type(hi), int)
                self.assertLessEqual(
                    lo, min(v[key] for v in vectors), (order, edges, key))
                self.assertLessEqual(
                    max(v[key] for v in vectors), hi, (order, edges, key))


class StateConeTests(unittest.TestCase):
    def test_state_intervals_contain_exact_catalog_stratum_pairs(self):
        # The three m=3 fixture stratum pairs; exact h and g values from
        # the naive reference vectors of real catalog graphs.
        # CLAMP: the first <= 20 catalog graphs per side (the order-24
        # 132-edge side holds exactly two; the order-23 122-edge side
        # holds exactly two), one cached 16 MB census scan in total.
        edge_bounds = load_edge_bounds(DATA / "structural_tables.json")
        maxima = {o: bounds[1] for o, bounds in edge_bounds.items()}
        pairs = (
            (20, 100, "r4520.100.g6", 24, 132, None),
            (21, 107, "r4521.107.g6", 23, 122, "r4523.122.g6"),
            (22, 114, "r4522.114.g6", 22, 114, "r4522.114.g6"),
        )
        for d, ex, x_name, m, ey, y_name in pairs:
            self.assertEqual(m, 44 - d)
            x_path = EXTREME / x_name
            y_path = DATA / "r45_24.g6" if y_name is None else EXTREME / y_name
            state = m4_state_interval(d, maxima[d] - ex, maxima[m] - ey)
            xs = read_graphs(x_path, ex, 20)
            ys = read_graphs(y_path, ey, 20)
            self.assertTrue(xs)
            self.assertTrue(ys)
            self.assertEqual(len(xs[0]), d)
            self.assertEqual(len(ys[0]), m)
            xv = [naive_motif4_vector(adj) for adj in xs]
            yv = [naive_motif4_vector(adj) for adj in ys]
            for vx, vy in itertools.product(xv, yv):
                exact_h = (3 * (47 - 2 * d) * vx["t"] + 4 * vx["diamond"]
                           - 6 * vx["k3k1"] - 12 * vy["i4"])
                self.assertLessEqual(state.h_lo, exact_h, (d, ex))
                self.assertLessEqual(exact_h, state.h_hi, (d, ex))
                exact_g = (3 * (c3(m) - ey * (m - 2) + vy["q"])
                           - (48 - 3 * d) * ex - 3 * vx["q"])
                self.assertLessEqual(state.g_lo, exact_g, (d, ex))
                self.assertLessEqual(exact_g, state.g_hi, (d, ex))

    def test_build_states_reproduces_frozen_m3_domain_and_records(self):
        # CLAMP: the frozen m=3 artifact is parsed once (module cache).
        states = build_states()
        self.assertEqual(len(states), 3215)
        self.assertEqual(states, sorted(states))
        records = m3_artifact_states()
        self.assertEqual(len(records), 3215)
        seven = ("d", "a", "b", "deficiency", "excess_balance",
                 "g_lo", "g_hi")
        for state, record in zip(states, records):
            self.assertIsInstance(state, State)
            for field in seven:
                self.assertEqual(getattr(state, field), record[field],
                                 (state.d, state.a, state.b, field))
            self.assertLessEqual(state.h_lo, state.h_hi)
        self.assertEqual({s.d for s in states}, set(STATE_DEGREES))
        self.assertEqual(
            tuple(f.name for f in dataclasses.fields(State)), STATE_FIELDS)
        with self.assertRaises(dataclasses.FrozenInstanceError):
            states[0].h_lo = 0   # the dataclass is frozen

    def test_state_record_serializes_the_v2_schema(self):
        states = build_states()
        for state in (states[0], states[-1]):
            record = state_record(state)
            self.assertEqual(tuple(sorted(record)),
                             tuple(sorted(STATE_RECORD_KEYS)))
            for field in STATE_FIELDS:
                self.assertIs(type(record[field]), int)
                self.assertEqual(record[field], getattr(state, field))
            for key in ("x_interval_source", "y_interval_source",
                        "x_motif_source", "y_motif_source"):
                self.assertIsInstance(record[key], str)
                self.assertTrue(record[key])


class KernelReplayTests(unittest.TestCase):
    def test_kernel_replays_656_residuals_and_n49_ground_truth(self):
        # CLAMP: all 328 published R(5,5,42) graphs once (656 residual
        # evaluations through the landed Task 1 kernel, ~10-15 s); the
        # n=49 replay reuses the cached two 132-edge census graphs.
        published = read_graphs(DATA / "r55_42some.g6")
        self.assertEqual(len(published), 328)
        tested = 0
        for adj in published:
            self.assertEqual(m4_residual_specialized(adj), 0)
            self.assertEqual(m4_residual_specialized(complement(adj)), 0)
            tested += 2
        self.assertEqual(tested, 656)

        extremal = read_graphs(DATA / "r45_24.g6", 132, 20)
        self.assertEqual(len(extremal), 2)
        for x in extremal:
            counts = motif4(x)
            self.assertEqual(counts.triangles, 176)
            self.assertEqual(counts.diamonds, 792)
            self.assertEqual(counts.k3k1, 528)
        self.assertEqual(
            sorted(independent_quad_count(x) for x in extremal),
            [138, 144],
        )
        row = 3 * (49 + 2 - 2 * 24) * 176 + 4 * 792 - 6 * 528
        self.assertEqual(row, 1584)      # = 12 * 132, the forced i4 mean
        self.assertEqual(row, 12 * 132)
        for x in extremal:
            self.assertGreater(independent_quad_count(x), 132)


class RelaxedEnvelopeTests(unittest.TestCase):
    def test_relaxed_rhs_intervals_equal_naive_histogram_aggregates(self):
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        order = 6   # CLAMP: every edge class of one small order.
        compared = 0
        for edges in range(c2(order) + 1):
            expected = naive_relaxed_rhs(order, edges, windows)
            if expected is None:
                continue
            self.assertEqual(
                relaxed_rhs_intervals(order, edges), expected,
                (order, edges),
            )
            compared += 1
        self.assertGreaterEqual(compared, 10)
        # R2's right-hand side is class-determined and stays exact.
        self.assertEqual(relaxed_rhs_intervals(6, 7)["R2"], (42, 42))

    def test_verified_basic_solution_count_is_independently_exact(self):
        # CLAMP: one pinned small class, (5, 5); the in-test fraction-free
        # enumeration of all C(17, 9) = 24,310 candidate bases is the
        # suite's largest single cost (~5-15 s).  No wall time recorded.
        order, edges = 5, 5
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        rhs = relaxed_rhs_intervals(order, edges)
        self.assertEqual(rhs, naive_relaxed_rhs(order, edges, windows))
        count, exact = enumerate_basic_solutions_exact(order, edges, rhs)
        self.assertGreater(count, 0)
        self.assertLessEqual(count, 24310)
        self.assertEqual(verified_basic_solution_count(order, edges), count)
        for key in MOTIF4_KEYS:
            lo, hi = exact[key]
            self.assertEqual(
                outer_motif_window(order, edges, key),
                (math.floor(lo), math.ceil(hi)), key,
            )
        # C5 is a real R(4,5,5) graph of the pinned class; every one of
        # its aggregate values must lie inside the exact windows.
        cycle5 = [18, 5, 10, 20, 9]
        vector = naive_motif4_vector(cycle5)
        self.assertEqual(vector,
                         {"q": 5, "t": 0, "diamond": 0, "k3k1": 0, "i4": 0})
        for key in MOTIF4_KEYS:
            self.assertLessEqual(exact[key][0], vector[key], key)
            self.assertLessEqual(vector[key], exact[key][1], key)

    def test_enumeration_failure_falls_back_to_sound_trivial_window(self):
        # CLAMP: one class, (6, 13), distinct from every other class the
        # suite queries, so no internal LP cache can mask the fallback.
        order, edges = 6, 13
        with mock.patch.object(
                m4, "enumerate_basic_solutions",
                side_effect=RuntimeError("forced certification failure")):
            for key in MOTIF4_KEYS:
                self.assertEqual(outer_motif_window(order, edges, key),
                                 (0, c4(order)), key)
        self.assertIn((order, edges), m4.lp_fallback_classes())


if __name__ == "__main__":
    unittest.main()
