#!/usr/bin/env python3
"""Exact mixed-identity motif windows + extended n=45 state cone (Task 2 tests).

Every observable contract below is specified in
notes/engstrom_identity_mixed_plan_2026-08-18.md, Task 2.
"""

import itertools
import json
import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
sys.path.insert(0, str(ROOT / "src"))

from check_ramsey import parse_graph6_line  # noqa: E402
from subgraph_identities import complement, induced_subgraph  # noqa: E402
from mixed_subgraph_identities import (  # noqa: E402
    mixed_residual,
    mixed_vertex_row,
    motif_sweep,
)
from m4_deficiency_cone import (  # noqa: E402 (frozen, re-imported)
    affine_interval,
    enumerate_basic_solutions,
)
from mixed_deficiency_cone import (  # noqa: E402  (RED before Task 2)
    CAMPAIGN_ID,
    EXPECTED_STATES,
    MOTIF4_KEYS,
    SCHEMA_VERSION,
    STREAM_KEYS,
    State,
    _enumerate_basic_solutions,
    build_states,
    catalog_motif_windows,
    lp_fallback_classes,
    mixed_state_interval,
    motif_interval,
    outer_motif_window,
    state_record,
    stream_catalog_motifs,
)


# ---------------------------------------------------- naive reference -----

_PB = {(i, j): bit for bit, (i, j)
       in enumerate(itertools.combinations(range(4), 2))}


def _canon(tuple_rows):
    best = None
    for perm in itertools.permutations(range(4)):
        pat = 0
        for i in range(4):
            for j in range(i + 1, 4):
                if (tuple_rows[perm[i]] >> perm[j]) & 1:
                    pat |= 1 << _PB[(i, j)]
        if best is None or pat < best:
            best = pat
    return best


_CANONPAT = {}
for _p in range(1 << 6):
    _rows = [0] * 4
    for (u, v), bit in _PB.items():
        if (_p >> bit) & 1:
            _rows[u] |= 1 << v
            _rows[v] |= 1 << u
    _CANONPAT[_p] = _canon(tuple(_rows))


def _class_pattern(edges):
    a = [0] * 4
    for (i, j) in edges:
        a[i] |= 1 << j
        a[j] |= 1 << i
    return _canon(tuple(a))


_CLS = {
    "i4": _class_pattern([]),
    "k2_2k1": _class_pattern([(0, 1)]),
    "two_k2": _class_pattern([(0, 1), (2, 3)]),
    "p3_k1": _class_pattern([(0, 1), (1, 2)]),
    "p4": _class_pattern([(0, 1), (1, 2), (2, 3)]),
    "claw": _class_pattern([(0, 1), (0, 2), (0, 3)]),
    "k3k1": _class_pattern([(0, 1), (0, 2), (1, 2)]),
    "c4": _class_pattern([(0, 1), (0, 2), (1, 3), (2, 3)]),
    "paw": _class_pattern([(0, 1), (0, 2), (0, 3), (1, 2)]),
    "diamond": _class_pattern([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]),
    "k4": _class_pattern([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
}


def naive_sweep_vector(adj):
    """The full 15-key aggregate vector by 4-subset + 3-subset enumeration."""
    n = len(adj)
    raw = {key: 0 for key in STREAM_KEYS}
    wedges = sum(r.bit_count() * (r.bit_count() - 1) // 2 for r in adj)
    edges = sum(r.bit_count() for r in adj) // 2
    tpairs = 0
    for a in range(n):
        for b in range(a + 1, n):
            if not ((adj[a] >> b) & 1):
                continue
            tpairs += (adj[a] & adj[b]).bit_count()
    t = tpairs // 3
    raw["t"] = t
    raw["q"] = wedges - t
    for quad in itertools.combinations(range(n), 4):
        pat = 0
        for (i, j), bit in _PB.items():
            if (adj[quad[i]] >> quad[j]) & 1:
                pat |= 1 << bit
        cls = _CANONPAT[pat]
        for key, cp in _CLS.items():
            if cls == cp:
                raw[key] += 1
                break
    raw["p3"] = wedges - 3 * t
    raw["pc3"] = edges * (n - 2) - 2 * raw["p3"] - 3 * t
    return raw


def _c2(value):
    return value * (value - 1) // 2


def _c3(value):
    return value * (value - 1) * (value - 2) // 6


def sweep_vector(adj):
    sw = motif_sweep(adj)
    return {
        "q": sw.wedges - sw.triangles,
        "t": sw.triangles,
        "p3": sw.induced_p3,
        "pc3": sw.pc3,
        "i4": sw.e4,
        "k4": sw.k4,
        "diamond": sw.diamond,
        "k3k1": sw.k3k1,
        "c4": sw.c4,
        "claw": sw.claw,
        "paw": sw.paw,
        "p4": sw.p4,
        "two_k2": sw.two_k2,
        "k2_2k1": sw.k2_2k1,
        "p3_k1": sw.p3_k1,
    }


def _edges_of(adj):
    return sum(r.bit_count() for r in adj) // 2


FIXTURE_FILES = {
    (17, 77): DATA / "r45extreme" / "r4517.77.g6",
    (20, 100): DATA / "r45extreme" / "r4520.100.g6",
    (21, 107): DATA / "r45extreme" / "r4521.107.g6",
    (22, 114): DATA / "r45extreme" / "r4522.114.g6",
    (23, 122): DATA / "r45extreme" / "r4523.122.g6",
}


def _class_graphs(cls):
    order, edges = cls
    out = []
    with FIXTURE_FILES[cls].open(encoding="ascii") as stream:
        for line in stream:
            if not line.strip():
                continue
            o, adj = parse_graph6_line(line)
            assert o == order and _edges_of(adj) == edges
            out.append(adj)
    return out


def _order24_132_graphs():
    out = []
    with (DATA / "r45_24.g6").open(encoding="ascii") as stream:
        for line in stream:
            text = line.strip()
            if not text:
                continue
            order, adj = parse_graph6_line(text)
            if _edges_of(adj) == 132:
                out.append(adj)
    return out


class MixedConeTests(unittest.TestCase):
    def test_sweep_vector_matches_naive_on_catalog_slice(self):
        graphs = []
        with FIXTURE_FILES[(17, 77)].open(encoding="ascii") as stream:
            for line in stream:
                if not line.strip():
                    continue
                _, adj = parse_graph6_line(line)
                graphs.append(adj)
            self.assertEqual(len(graphs), 7060)
        for adj in graphs[:250]:
            got = sweep_vector(adj)
            want = naive_sweep_vector(adj)
            self.assertEqual(set(got), set(STREAM_KEYS))
            self.assertEqual(set(want), set(STREAM_KEYS))
            self.assertEqual(got, want)
        for adj in _order24_132_graphs():
            self.assertEqual(sweep_vector(adj), naive_sweep_vector(adj))

    def test_affine_interval_repinned(self):
        self.assertEqual(affine_interval(-3, (10, 20)), (-60, -30))
        self.assertEqual(affine_interval(4, (-5, 2)), (-20, 8))
    def test_path_only_stream_contract(self):
        records = stream_catalog_motifs(FIXTURE_FILES[(17, 77)])
        self.assertEqual(set(records), {(17, 77)})
        self.assertEqual(set(records[(17, 77)]), set(STREAM_KEYS))
        self.assertEqual(sum(records[(17, 77)]["q"].values()), 7060)
    def test_mixed_lp_projection_matches_m4(self):
        m4_count, m4_best = enumerate_basic_solutions(6, 13)
        mixed_count, mixed_best = _enumerate_basic_solutions(6, 13)
        self.assertEqual(mixed_count, m4_count)
        for key in MOTIF4_KEYS:
            self.assertEqual(mixed_best[key], m4_best[key], key)



    def test_stream_and_windows_override_outer(self):
        path = FIXTURE_FILES[(22, 114)]
        graphs, records = stream_catalog_motifs(path, 22, 114, 133)
        self.assertEqual(graphs, 133)
        windows = catalog_motif_windows(records)
        self.assertEqual(set(windows[(22, 114)]), set(STREAM_KEYS))
        catalog_map = {cls: windows[cls] for cls in windows}
        for key in STREAM_KEYS:
            self.assertEqual(
                motif_interval(22, 114, key, catalog_map),
                windows[(22, 114)][key],
            )
            # exact catalog windows override the (sound, possibly wider) outer
            window = outer_motif_window(22, 114, key)
            lo, hi = window
            obs = windows[(22, 114)][key]
            self.assertLessEqual(lo, obs[0], key)
            self.assertLessEqual(obs[1], hi, key)
            # outer need not collapse to the catalog window (override check)
        # a key absent from the map resolves to the outer window
        self.assertEqual(
            motif_interval(22, 114, "t", {}),
            outer_motif_window(22, 114, "t"),
        )

    def test_outer_windows_cover_catalog_values(self):
        for cls, path in FIXTURE_FILES.items():
            order, edges = cls
            graphs, records = stream_catalog_motifs(path, order, edges, None)
            windows = catalog_motif_windows(records)[(order, edges)]
            for key in STREAM_KEYS:
                lo, hi = outer_motif_window(order, edges, key)
                self.assertLessEqual(lo, windows[key][0], (cls, key))
                self.assertLessEqual(windows[key][1], hi, (cls, key))

    def test_fixture_state_intervals_contain_exact_rows(self):
        catalog = {}
        for cls in ((22, 114),):
            _, records = stream_catalog_motifs(
                FIXTURE_FILES[cls], cls[0], cls[1], None)
            catalog[cls] = catalog_motif_windows(records)[cls]
        _, records = stream_catalog_motifs(
            FIXTURE_FILES[(20, 100)], 20, 100, None)
        catalog[(20, 100)] = catalog_motif_windows(records)[(20, 100)]
        _, records = stream_catalog_motifs(
            FIXTURE_FILES[(21, 107)], 21, 107, None)
        catalog[(21, 107)] = catalog_motif_windows(records)[(21, 107)]
        _, records = stream_catalog_motifs(
            FIXTURE_FILES[(23, 122)], 23, 122, None)
        catalog[(23, 122)] = catalog_motif_windows(records)[(23, 122)]
        z24 = _order24_132_graphs()
        z24_records = {(24, 132): {key: {} for key in STREAM_KEYS}}
        for adj in z24:
            vec = sweep_vector(adj)
            for key in STREAM_KEYS:
                level = z24_records[(24, 132)][key]
                level[vec[key]] = level.get(vec[key], 0) + 1
        catalog[(24, 132)] = catalog_motif_windows(z24_records)[(24, 132)]
        x_graphs = {
            (20, 100): _class_graphs((20, 100)),
            (21, 107): _class_graphs((21, 107)),
            (22, 114): _class_graphs((22, 114)),
        }
        z_graphs = {
            (24, 132): z24,
            (23, 122): _class_graphs((23, 122)),
            (22, 114): _class_graphs((22, 114)),
        }
        e45 = {20: 100, 21: 107, 22: 114, 23: 122, 24: 132}
        for (x_cls, z_cls) in (
            ((20, 100), (24, 132)),
            ((21, 107), (23, 122)),
            ((22, 114), (22, 114)),
        ):
            d = x_cls[0]
            m = z_cls[0]
            a = e45[d] - x_cls[1]
            b = e45[m] - z_cls[1]
            state = mixed_state_interval(d, a, b, catalog, {})
            self.assertLessEqual(state.f_lo, state.f_hi)
            for x in x_graphs[x_cls][:20]:
                for z in z_graphs[z_cls][:20]:
                    exact = mixed_vertex_row(
                        45, d, motif_sweep(x), motif_sweep(z))
                    self.assertLessEqual(state.f_lo, exact, (x_cls, z_cls))
                    self.assertLessEqual(exact, state.f_hi, (x_cls, z_cls))
            self.assertLessEqual(state.g_lo, state.g_hi)
            self.assertLessEqual(state.h_lo, state.h_hi)

    def test_cone_stability_and_frozen_projection(self):
        artifact = json.loads(
            (DATA / "higher_identity_m4.json").read_text(encoding="ascii")
        )
        frozen = {
            (s["d"], s["a"], s["b"]): s for s in artifact["states"]
        }
        states = build_states()
        self.assertEqual(len(states), EXPECTED_STATES)
        seen = {}
        for state in states:
            self.assertLessEqual(state.f_lo, state.f_hi)
            self.assertLessEqual(state.g_lo, state.g_hi)
            self.assertLessEqual(state.h_lo, state.h_hi)
            key = (state.d, state.a, state.b)
            self.assertIn(key, frozen)
            self.assertNotIn(key, seen)
            seen[key] = state
            rec = frozen[key]
            self.assertEqual(state.g_lo, rec["g_lo"], key)
            self.assertEqual(state.g_hi, rec["g_hi"], key)
            self.assertEqual(state.h_lo, rec["h_lo"], key)
            self.assertEqual(state.h_hi, rec["h_hi"], key)
            self.assertEqual(state.deficiency, rec["deficiency"], key)
            self.assertEqual(state.excess_balance, rec["excess_balance"], key)
            rebuilt = state_record(state)
            self.assertEqual(
                rebuilt["x_interval_source"], rec["x_interval_source"], key
            )
            self.assertEqual(
                rebuilt["y_interval_source"], rec["y_interval_source"], key
            )

    def test_relation_rows_hold_on_random_graphs(self):
        # Exact per-vertex row-space identities (the frozen m4 wiring).
        # R1: sum_i s_i = C(n,4); R2: sum_i e_i s_i = e*C(n-2,2).
        # R3: sum_i e_i^2 s_i = e*C(n-2,2) + 2W(n-4) + e(e-1), W = wedges.
        # R4: sum_i sigma_i s_i = sum_v [d*C(n-1-d,2)+4*C(d,2)(n-1-d)+9*C(d,3)].
        # R5: sum_i tau_i s_i = (n-3)*t.
        from m4_deficiency_cone import relation_rows
        rows = relation_rows()
        rng = random.Random(20260818)
        for _ in range(60):
            n = rng.randint(2, 7)
            adj = [0] * n
            for u, v in itertools.combinations(range(n), 2):
                if rng.random() < 0.5:
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
            sw = motif_sweep(adj)
            order_map = ["e4", "k2_2k1", "two_k2", "p3_k1", "p4", "claw",
                         "k3k1", "c4", "paw", "diamond", "k4"]
            values = [getattr(sw, key) for key in order_map]
            deg = [r.bit_count() for r in adj]
            n = len(adj)
            n4 = n * (n - 1) * (n - 2) * (n - 3) // 24
            e = sw.edges
            wedges = sum(d * (d - 1) // 2 for d in deg)
            self.assertEqual(sum(values), n4, "R1")
            self.assertEqual(sum(r * v for r, v in zip(rows["R2"], values)),
                             e * _c2(n - 2), "R2")
            self.assertEqual(
                sum(r * v for r, v in zip(rows["R3"], values)),
                e * _c2(n - 2) + 2 * wedges * (n - 4) + e * (e - 1), "R3")
            self.assertEqual(
                sum(r * v for r, v in zip(rows["R4"], values)),
                sum(
                    d * _c2(n - 1 - d) + 4 * _c2(d) * (n - 1 - d) + 9 * _c3(d)
                    for d in deg
                ),
                "R4")
            self.assertEqual(
                sum(r * v for r, v in zip(rows["R5"], values)),
                (n - 3) * sw.triangles, "R5")

    def test_replays(self):
        published = []
        with (ROOT / "data" / "r55_42some.g6").open(encoding="ascii") as fh:
            for line in fh:
                if line.strip():
                    _, adj = parse_graph6_line(line)
                    published.append(adj)
        self.assertEqual(len(published), 328)
        for adj in published:
            self.assertEqual(mixed_residual(adj), 0)
            self.assertEqual(mixed_residual(complement(adj)), 0)
        z24 = _order24_132_graphs()
        self.assertEqual(len(z24), 2)
        for z in z24:
            sw = motif_sweep(z)
            self.assertEqual(sw.triangles, 176)
            self.assertEqual(sw.diamond, 792)
            self.assertEqual(sw.paw, 1584)
            self.assertEqual(sw.claw, 792)
            self.assertEqual(sw.k3k1, 528)
        values = sorted(
            mixed_vertex_row(49, 24, motif_sweep(xi), motif_sweep(zj))
            for xi in z24 for zj in z24
        )
        self.assertEqual(values, [0, 144, 288, 432])
        self.assertEqual(sum(v == 0 for v in values), 1)

    def test_campaign_constants(self):
        self.assertEqual(SCHEMA_VERSION, 3)
        self.assertEqual(CAMPAIGN_ID,
                         "higher_order_identity_positive_deficiency_mixed")
        self.assertEqual(EXPECTED_STATES, 3215)


if __name__ == "__main__":
    unittest.main()
