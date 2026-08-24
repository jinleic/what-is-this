#!/usr/bin/env python3
"""Exact counterexamples to the local charging step in Clement--Bader (2007).

The draft claims that a point incident with more than two lines belongs to at
most two pairs of triangular faces sharing a side.  Both examples below use
only triangular arrangement faces.  The triple flower has six shared sides at
one triple point; the quadruple flower has eight at one quadruple point and
saturates the 2k local-token term in the campaign's capacity theorem.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "math" / "kobon"))
import engine  # noqa: E402


def shared_sides(lines, selected):
    ms = [F(m) for m, _ in lines]
    bs = [F(b) for _, b in lines]
    X = engine.crossings(len(lines), ms, bs)
    uses = Counter()
    for t in selected:
        for line in t:
            others = [x for x in t if x != line]
            endpoints = tuple(sorted((
                engine.cross(X, line, others[0]),
                engine.cross(X, line, others[1]),
            )))
            uses[line, endpoints] += 1
    return uses


def verify(name, lines, selected, multiplicity, expected_shared):
    n = len(lines)
    ms = [F(m) for m, _ in lines]
    bs = [F(b) for _, b in lines]
    selected = [tuple(t) for t in selected]
    ok, reason = engine.verify_selection(n, ms, bs, selected, len(selected))
    if not ok:
        raise AssertionError(f"{name}: selection failed: {reason}")
    crossings, bounded_faces = engine.face_multiplicity_counts(
        n, ms, bs, selected)
    if crossings != 0:
        raise AssertionError(f"{name}: selected triangle is not a face")

    point_lines = {}
    X = engine.crossings(n, ms, bs)
    for i, j in combinations(range(n), 2):
        p = engine.cross(X, i, j)
        if p is not None:
            point_lines.setdefault(p, set()).update((i, j))
    origin_lines = point_lines.get((F(0), F(0)), set())
    if len(origin_lines) != multiplicity:
        raise AssertionError(
            f"{name}: origin multiplicity {len(origin_lines)} != {multiplicity}")

    uses = shared_sides(lines, selected)
    shared = [(line, endpoints) for (line, endpoints), count in uses.items()
              if count == 2]
    if len(shared) != expected_shared or max(uses.values()) != 2:
        raise AssertionError(
            f"{name}: shared={len(shared)}, max_use={max(uses.values())}")
    if any(line not in origin_lines or (F(0), F(0)) not in endpoints
           for line, endpoints in shared):
        raise AssertionError(f"{name}: a shared side is not radial at the origin")

    q = sum(ms[i] == ms[j] for i, j in combinations(range(n), 2))
    return {
        "name": name,
        "n": n,
        "lines_y_equals_mx_plus_b": [
            [str(m), str(b)] for m, b in zip(ms, bs)
        ],
        "selected_triangular_faces": [list(t) for t in selected],
        "selected_count": len(selected),
        "line_triangle_crossing_incidences": crossings,
        "bounded_arrangement_faces": bounded_faces,
        "parallel_pairs": q,
        "origin_multiplicity": len(origin_lines),
        "shared_elementary_sides_at_origin": len(shared),
        "expected_local_token_maximum_2k": 2 * multiplicity,
        "capacity_point_coefficient_k_times_k_minus_4": (
            multiplicity * (multiplicity - 4)
        ),
        "verified": True,
    }


def main():
    triple = verify(
        "triple_hexagon_flower",
        [(-1, 0), (0, 0), (1, 0),
         (1, -2), (-1, 2), (0, 1), (1, 2), (-1, -2), (0, -1)],
        [(0, 1, 3), (1, 2, 4), (0, 2, 5),
         (0, 1, 6), (1, 2, 7), (0, 2, 8)],
        multiplicity=3,
        expected_shared=6,
    )
    quadruple = verify(
        "quadruple_triangle_flower",
        [(-1, 0), (0, 0), (1, 0), (2, 0),
         (-4, 6), (F(4, 3), F(2, 3)), (0, -2)],
        [(0, 1, 4), (0, 1, 5), (0, 3, 5), (0, 3, 6),
         (1, 2, 4), (1, 2, 5), (2, 3, 4), (2, 3, 6)],
        multiplicity=4,
        expected_shared=8,
    )
    print(json.dumps({"examples": [triple, quadruple]}, indent=2))


if __name__ == "__main__":
    main()
