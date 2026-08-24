#!/usr/bin/env python3
"""Exact rational counterexample to the face-specific square penalty.

The source arrangement is the classical simplicial arrangement A(12,1): the
six side lines and six symmetry axes of a regular hexagon.  We use a rational
affine image of the hexagon and then an exact generic projective chart, so all
12 resulting affine lines have finite, pairwise distinct slopes.

The replay uses two independent exact triangular-face predicates and checks the
multipoint census, bounded-face/edge identities, and shared-edge accounting.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import cmp_to_key
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
ENGINE_DIR = ROOT / "math" / "kobon"
if not (ENGINE_DIR / "engine.py").is_file():
    ENGINE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ENGINE_DIR))
import engine  # noqa: E402


HEXAGON = (
    (F(2), F(0)), (F(1), F(1)), (F(-1), F(1)),
    (F(-2), F(0)), (F(-1), F(-1)), (F(1), F(-1)),
)


def line_through(p, q):
    x1, y1 = p
    x2, y2 = q
    return y1 - y2, x2 - x1, x1 * y2 - x2 * y1


def midpoint(p, q):
    return (p[0] + q[0]) / 2, (p[1] + q[1]) / 2


def projective_a12_1_lines():
    """Return y = m*x+b lines in a generic exact chart of A(12,1)."""
    coefficients = [
        line_through(HEXAGON[i], HEXAGON[(i + 1) % 6])
        for i in range(6)
    ]
    coefficients.extend(
        line_through(HEXAGON[i], HEXAGON[i + 3])
        for i in range(3)
    )
    coefficients.extend(
        line_through(
            midpoint(HEXAGON[i], HEXAGON[(i + 1) % 6]),
            midpoint(HEXAGON[i + 3], HEXAGON[(i + 4) % 6]),
        )
        for i in range(3)
    )

    # First apply (X,Y)=(x+2y,3x+5y), then choose the new line at infinity
    # Z + alpha*X + beta*Y = 0.  For A0*X+B0*Y+C*Z=0 this gives
    # (A0-alpha*C)*X + (B0-beta*C)*Y + C*Z_new = 0.
    alpha, beta = F(1, 17), F(7, 19)
    lines = []
    for a, b, c in coefficients:
        a0, b0 = -5 * a + 3 * b, 2 * a - b
        a1, b1 = a0 - alpha * c, b0 - beta * c
        if b1 == 0:
            raise AssertionError("projective chart produced a vertical line")
        lines.append((-a1 / b1, -c / b1))
    if len(lines) != 12 or len(set(lines)) != 12:
        raise AssertionError("construction did not produce 12 distinct lines")
    return tuple(lines)


def intersection(line1, line2):
    m1, b1 = line1
    m2, b2 = line2
    if m1 == m2:
        return None
    x = (b2 - b1) / (m1 - m2)
    return x, m1 * x + b1


def triangle_vertices(lines, triple):
    i, j, k = triple
    vertices = (
        intersection(lines[i], lines[j]),
        intersection(lines[i], lines[k]),
        intersection(lines[j], lines[k]),
    )
    if None in vertices or len(set(vertices)) != 3:
        return None
    return vertices


def straddle_face(lines, triple):
    """A triangle is a face iff no outside line crosses its open interior."""
    vertices = triangle_vertices(lines, triple)
    if vertices is None:
        return False
    for outside, (m, b) in enumerate(lines):
        if outside in triple:
            continue
        signs = [y - m * x - b for x, y in vertices]
        if min(signs) < 0 < max(signs):
            return False
    return True


def direct_gap_face(lines, triple):
    """Independently reject any strict outside intersection on a side."""
    if triangle_vertices(lines, triple) is None:
        return False
    for side in triple:
        others = [line for line in triple if line != side]
        endpoints = [intersection(lines[side], lines[line]) for line in others]
        x0, x1 = sorted((endpoints[0][0], endpoints[1][0]))
        for outside in range(len(lines)):
            if outside in triple:
                continue
            point = intersection(lines[side], lines[outside])
            if point is not None and x0 < point[0] < x1:
                return False
    return True

def polar_ray_order(lines, incident):
    """Sort the two directed rays of each incident line by exact polar angle."""
    rays = [(line, sign) for line in incident for sign in (1, -1)]

    def vector(ray):
        line, sign = ray
        return F(sign), F(sign) * lines[line][0]

    def half(vector):
        x, y = vector
        return 0 if y > 0 or (y == 0 and x > 0) else 1

    def compare(left, right):
        u, v = vector(left), vector(right)
        hu, hv = half(u), half(v)
        if hu != hv:
            return -1 if hu < hv else 1
        cross = u[0] * v[1] - u[1] * v[0]
        if cross:
            return -1 if cross > 0 else 1
        raise AssertionError("two incident lines produced the same ray")

    return sorted(rays, key=cmp_to_key(compare))


def multipoint_sector_census(lines, faces, point_lines, edge_uses):
    """Compute the exact triangular-sector run identity at every multipoint."""
    records = []
    for point, incident in sorted(point_lines.items()):
        k = len(incident)
        if k < 3:
            continue
        rays = polar_ray_order(lines, incident)
        ray_index = {ray: index for index, ray in enumerate(rays)}
        sectors = [False] * (2 * k)

        for triple in faces:
            through = [line for line in triple if line in incident]
            if len(through) != 2:
                continue
            third = next(line for line in triple if line not in through)
            face_rays = []
            for side in through:
                other_vertex = intersection(lines[side], lines[third])
                if other_vertex == point:
                    raise AssertionError("nondegenerate face collapsed at a multipoint")
                sign = 1 if other_vertex[0] > point[0] else -1
                face_rays.append((side, sign))
            left, right = (ray_index[ray] for ray in face_rays)
            if (left + 1) % len(rays) == right:
                sector = left
            elif (right + 1) % len(rays) == left:
                sector = right
            else:
                raise AssertionError("triangular face does not occupy one local sector")
            if sectors[sector]:
                raise AssertionError("two faces occupy the same local sector")
            sectors[sector] = True

        triangular_sectors = sum(sectors)
        nontriangular_sectors = len(sectors) - triangular_sectors
        runs = (
            0 if nontriangular_sectors == 0
            else sum(value and not sectors[index - 1]
                     for index, value in enumerate(sectors))
        )
        shared_incident = sum(
            use == 2 and point in endpoints
            for (_, endpoints), use in edge_uses.items()
        )
        predicted_shared_incident = 2 * k - nontriangular_sectors - runs
        if shared_incident != predicted_shared_incident:
            raise AssertionError("triangular-sector run identity failed")
        local_defect = shared_incident - 2 * (k - 2)
        if local_defect != 4 - nontriangular_sectors - runs:
            raise AssertionError("local defect identity failed")
        records.append({
            "point": [str(point[0]), str(point[1])],
            "multiplicity": k,
            "triangular_sectors": triangular_sectors,
            "nontriangular_sectors": nontriangular_sectors,
            "triangular_sector_runs": runs,
            "incident_shared_edges": shared_incident,
            "local_defect_4_minus_z_minus_r": local_defect,
        })
    return records


def verify():
    lines = projective_a12_1_lines()
    n = len(lines)
    triples = list(combinations(range(n), 3))
    straddle_faces = {triple for triple in triples if straddle_face(lines, triple)}
    gap_faces = {triple for triple in triples if direct_gap_face(lines, triple)}
    if straddle_faces != gap_faces:
        raise AssertionError("the independent face predicates disagree")

    point_lines = {}
    line_points = {line: set() for line in range(n)}
    parallel_pairs = 0
    for i, j in combinations(range(n), 2):
        point = intersection(lines[i], lines[j])
        if point is None:
            parallel_pairs += 1
            continue
        point_lines.setdefault(point, set()).update((i, j))
        line_points[i].add(point)
        line_points[j].add(point)
    multiplicity_counts = Counter(len(incident) for incident in point_lines.values())
    multipoint_counts = {
        k: count for k, count in sorted(multiplicity_counts.items()) if k >= 3
    }

    edge_uses = Counter()
    for triple in gap_faces:
        for side in triple:
            others = [line for line in triple if line != side]
            endpoints = tuple(sorted((
                intersection(lines[side], lines[others[0]]),
                intersection(lines[side], lines[others[1]]),
            )))
            edge_uses[side, endpoints] += 1
    if any(use not in (1, 2) for use in edge_uses.values()):
        raise AssertionError("a face edge has impossible incidence")

    bounded_edges = sum(len(points) - 1 for points in line_points.values())
    used_once = sum(use == 1 for use in edge_uses.values())
    shared = sum(use == 2 for use in edge_uses.values())
    unused = bounded_edges - used_once - shared
    multipoints = {
        point for point, incident in point_lines.items() if len(incident) >= 3
    }
    shared_with_two_multipoint_endpoints = sum(
        use == 2 and endpoints[0] in multipoints and endpoints[1] in multipoints
        for (_, endpoints), use in edge_uses.items()
    )
    sector_census = multipoint_sector_census(
        lines, gap_faces, point_lines, edge_uses)
    sector_defect_sum = sum(
        record["local_defect_4_minus_z_minus_r"] for record in sector_census)
    if 3 * len(gap_faces) != used_once + 2 * shared:
        raise AssertionError("face-side incidence identity failed")

    square_penalty = sum(
        count * (k - 2) ** 2 for k, count in multipoint_counts.items()
    )
    linear_defect = sum(
        count * (k - 2) for k, count in multipoint_counts.items()
    )
    square_rhs = n * (n - 2) - 2 * parallel_pairs - square_penalty
    square_lhs = 3 * len(gap_faces)
    bounded_faces = (
        comb(n - 1, 2)
        - parallel_pairs
        - sum(count * comb(k - 1, 2)
              for k, count in multipoint_counts.items())
    )
    edge_formula = (
        n * (n - 2)
        - 2 * parallel_pairs
        - sum(count * k * (k - 2)
              for k, count in multipoint_counts.items())
    )

    if parallel_pairs != 0:
        raise AssertionError("generic chart retained a parallel pair")
    if multipoint_counts != {3: 15, 6: 1}:
        raise AssertionError(f"unexpected multipoint census {multipoint_counts}")
    if len(gap_faces) != 30 or bounded_faces != 30:
        raise AssertionError("A(12,1) did not have 30 triangular bounded faces")
    if bounded_edges != edge_formula or (used_once, shared, unused) != (12, 39, 0):
        raise AssertionError("bounded-edge census failed")
    if not (square_lhs == 90 and square_rhs == 89):
        raise AssertionError("counterexample margin changed")
    if shared - unused != 39 or 2 * linear_defect != 38:
        raise AssertionError("shared-minus-unused reformulation failed")
    if sum(record["incident_shared_edges"] for record in sector_census) != (
            shared + shared_with_two_multipoint_endpoints):
        raise AssertionError("shared-edge endpoint identity failed")
    if square_lhs - square_rhs != (
            sector_defect_sum - unused - shared_with_two_multipoint_endpoints):
        raise AssertionError("global sector-excess identity failed")

    ms, bs = zip(*lines)
    ok, reason = engine.verify_selection(
        n, list(ms), list(bs), sorted(gap_faces), minimum=30)
    if not ok:
        raise AssertionError(f"independent engine verification failed: {reason}")

    return {
        "status": "VERIFIED_COUNTEREXAMPLE",
        "source_arrangement": "A(12,1): regular-hexagon sides and symmetry axes",
        "regular_polygon_sides": 6,
        "symmetry_axis_types": {
            "opposite_vertex_axes": 3,
            "opposite_edge_midpoint_axes": 3,
        },
        "n": n,
        "lines_y_equals_mx_plus_b": [[str(m), str(b)] for m, b in lines],
        "parallel_pairs": parallel_pairs,
        "finite_vertex_multiplicities": {
            str(k): count for k, count in sorted(multiplicity_counts.items())
        },
        "multipoint_counts": {
            str(k): count for k, count in multipoint_counts.items()
        },
        "bounded_faces": bounded_faces,
        "triangular_faces": len(gap_faces),
        "triangular_face_triples": [list(triple) for triple in sorted(gap_faces)],
        "bounded_elementary_edges": bounded_edges,
        "face_edges_used_once": used_once,
        "face_edges_used_twice": shared,
        "unused_bounded_edges": unused,
        "linear_defect_sum": linear_defect,
        "shared_minus_unused": shared - unused,
        "twice_linear_defect": 2 * linear_defect,
        "shared_edges_with_two_multipoint_endpoints": (
            shared_with_two_multipoint_endpoints
        ),
        "multipoint_sector_census": sector_census,
        "sector_defect_sum": sector_defect_sum,
        "sector_excess_identity_rhs": (
            sector_defect_sum
            - unused
            - shared_with_two_multipoint_endpoints
        ),
        "square_penalty": square_penalty,
        "claimed_lhs_3F": square_lhs,
        "claimed_rhs": square_rhs,
        "violation_margin": square_lhs - square_rhs,
        "face_predicates": {
            "strict_straddle_count": len(straddle_faces),
            "direct_gap_count": len(gap_faces),
            "sets_equal": True,
        },
        "engine_verify_selection": reason,
    }


def main():
    result = verify()
    text = json.dumps(result, indent=2, sort_keys=True)
    if len(sys.argv) == 3 and sys.argv[1] == "--out":
        Path(sys.argv[2]).write_text(text + "\n")
    elif len(sys.argv) != 1:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} [--out FILE]")
    print(text)


if __name__ == "__main__":
    main()
