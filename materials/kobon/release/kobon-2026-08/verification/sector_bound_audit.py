#!/usr/bin/env python3
"""Exact-rational audit of face-sector and endpoint-closure bounds.

The audit samples deterministic arrangements with generic, parallel, and
concurrent degeneracies.  It checks the local facts encoded by
engine.add_selected_face_sector_bounds:

* at most four triangular faces use any fixed pair of lines;
* two faces extending along the same ray of a shared line have a common
  remote endpoint, hence force concurrency with their two third lines;
* the same endpoint theorem when two faces share only one supporting line
  and meet at a multipoint endpoint;
* two faces sharing an entire side segment lie on opposite sides;
* at a finite multipoint, at most two faces use a fixed pair;
* a pair separated on both arcs of the cyclic slope order supports no face.

Finite sampling checks the implementation and label-order convention; an
exact six-line control attains all four sectors simultaneously.
"""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
ENGINE_DIR = ROOT / "math" / "kobon"
SCRIPT_DIR = ROOT / "scratch" / "kobon"
if not (ENGINE_DIR / "engine.py").is_file():
    bundle = Path(__file__).resolve().parents[1]
    ENGINE_DIR = bundle / "scripts"
    SCRIPT_DIR = bundle / "scripts"
sys.path.insert(0, str(ENGINE_DIR))
sys.path.insert(0, str(SCRIPT_DIR))
import engine  # noqa: E402
from square_penalty_counterexample import direct_gap_face, intersection  # noqa: E402

MODES = (
    ("generic", 0, 0),
    ("parallel_pair", 1, 0),
    ("two_parallel_pairs", 2, 0),
    ("triple_point", 0, 1),
    ("two_triple_points", 0, 2),
    ("parallel_and_triple", 1, 1),
)


def sharpness_control() -> dict:
    lines = (
        (engine.Fr(0), engine.Fr(0)),
        (engine.Fr(1), engine.Fr(0)),
        (engine.Fr(2), engine.Fr(-2)),
        (engine.Fr(2, 3), engine.Fr(-2, 3)),
        (engine.Fr(2, 3), engine.Fr(2, 3)),
        (engine.Fr(2), engine.Fr(2)),
    )
    faces = sorted(
        triple
        for triple in combinations(range(6), 3)
        if direct_gap_face(lines, triple)
    )
    expected = [(0, 1, third) for third in range(2, 6)]
    if faces != expected:
        raise AssertionError(("sharpness faces", faces))

    vertex = intersection(lines[0], lines[1])
    keys = {}
    for _, _, third in faces:
        key = (
            intersection(lines[0], lines[third])[0] > vertex[0],
            intersection(lines[1], lines[third])[0] > vertex[0],
        )
        if key in keys:
            raise AssertionError(("sharpness sector collision", key))
        keys[key] = third
    if len(keys) != 4:
        raise AssertionError(("sharpness sector count", keys))
    return {
        "status": "PASS",
        "line_pair": [0, 1],
        "faces": [list(face) for face in faces],
        "sector_keys": {
            str(third): [int(first), int(second)]
            for (first, second), third in sorted(keys.items())
        },
    }


def audit(samples_per_mode: int) -> dict:
    counts = {
        "arrangements": 0,
        "triangular_faces": 0,
        "line_pair_checks": 0,
        "face_sector_assignments": 0,
        "shared_pair_face_pairs": 0,
        "shared_ray_cases": 0,
        "shared_pair_segment_face_pairs": 0,
        "multipoint_pair_checks": 0,
        "single_line_face_pairs": 0,
        "colliding_endpoint_cases": 0,
        "colliding_same_ray_cases": 0,
        "shared_segment_face_pairs": 0,
        "nonadjacent_multipoint_pair_checks": 0,
        "violations": 0,
    }
    examples = []
    for n in range(5, 11):
        for mode_index, (mode, parallel_pairs, triple_points) in enumerate(MODES):
            for sample in range(samples_per_mode):
                seed = 10_000_000 * n + 100_000 * mode_index + sample
                slopes, intercepts = engine.gen_config(
                    n,
                    n_par_pairs=parallel_pairs,
                    n_conc=triple_points,
                    coord=10**5,
                    seed=seed,
                )
                lines = tuple(zip(slopes, intercepts))
                faces = {
                    triple
                    for triple in combinations(range(n), 3)
                    if direct_gap_face(lines, triple)
                }
                point_lines = {}
                for i, j in combinations(range(n), 2):
                    point = intersection(lines[i], lines[j])
                    if point is not None:
                        point_lines.setdefault(point, set()).update((i, j))

                counts["arrangements"] += 1
                counts["triangular_faces"] += len(faces)
                for i, j in combinations(range(n), 2):
                    pair_faces = sorted(
                        face for face in faces if i in face and j in face
                    )
                    codegree = len(pair_faces)
                    counts["line_pair_checks"] += 1
                    if codegree > 4:
                        examples.append((mode, n, seed, i, j, "simple", codegree))

                    point = intersection(lines[i], lines[j])
                    if point is None:
                        continue

                    thirds = []
                    sector_keys = {}
                    for face in pair_faces:
                        third = next(line for line in face if line not in (i, j))
                        thirds.append(third)
                        key = (
                            intersection(lines[i], lines[third])[0] > point[0],
                            intersection(lines[j], lines[third])[0] > point[0],
                        )
                        counts["face_sector_assignments"] += 1
                        if key in sector_keys:
                            examples.append((
                                mode, n, seed, i, j, "sector_collision",
                                sector_keys[key], third,
                            ))
                        sector_keys[key] = third

                    for first, second in combinations(thirds, 2):
                        counts["shared_pair_face_pairs"] += 1
                        for shared in (i, j):
                            first_point = intersection(
                                lines[shared], lines[first])
                            second_point = intersection(
                                lines[shared], lines[second])
                            if ((first_point[0] > point[0])
                                    == (second_point[0] > point[0])):
                                counts["shared_ray_cases"] += 1
                                if first_point != second_point:
                                    examples.append((
                                        mode, n, seed, i, j, "shared_ray",
                                        shared, first, second,
                                    ))
                            if first_point == second_point:
                                counts[
                                    "shared_pair_segment_face_pairs"] += 1
                                pivot = j if shared == i else i
                                first_apex = intersection(
                                    lines[pivot], lines[first])
                                second_apex = intersection(
                                    lines[pivot], lines[second])
                                slope, intercept = lines[shared]
                                first_side = (
                                    first_apex[1]
                                    - slope * first_apex[0] - intercept)
                                second_side = (
                                    second_apex[1]
                                    - slope * second_apex[0] - intercept)
                                if first_side * second_side >= 0:
                                    examples.append((
                                        mode, n, seed, i, j,
                                        "shared_pair_segment_side",
                                        shared, first, second,
                                    ))

                    incident = point_lines[point]
                    if len(incident) < 3:
                        continue
                    counts["multipoint_pair_checks"] += 1
                    if codegree > 2:
                        examples.append((mode, n, seed, i, j, "multiple", codegree))

                    has_between = any(i < line < j for line in incident)
                    has_outside = any(line < i or line > j for line in incident)
                    if has_between and has_outside:
                        counts["nonadjacent_multipoint_pair_checks"] += 1
                        if codegree:
                            examples.append((
                                mode, n, seed, i, j, "nonadjacent", codegree))

                for shared in range(n):
                    shared_faces = sorted(
                        face for face in faces if shared in face
                    )
                    for first_face, second_face in combinations(
                            shared_faces, 2):
                        if set(first_face) & set(second_face) != {shared}:
                            continue
                        counts["single_line_face_pairs"] += 1
                        first_pair = tuple(
                            line for line in first_face if line != shared)
                        second_pair = tuple(
                            line for line in second_face if line != shared)
                        for first_endpoint in first_pair:
                            first_remote = next(
                                line for line in first_pair
                                if line != first_endpoint)
                            for second_endpoint in second_pair:
                                second_remote = next(
                                    line for line in second_pair
                                    if line != second_endpoint)
                                common_point = intersection(
                                    lines[shared], lines[first_endpoint])
                                if common_point != intersection(
                                        lines[shared], lines[second_endpoint]):
                                    continue
                                counts["colliding_endpoint_cases"] += 1
                                first_remote_point = intersection(
                                    lines[shared], lines[first_remote])
                                second_remote_point = intersection(
                                    lines[shared], lines[second_remote])
                                same_ray = (
                                    first_remote_point[0] > common_point[0]
                                ) == (
                                    second_remote_point[0] > common_point[0]
                                )
                                if not same_ray:
                                    continue
                                counts["colliding_same_ray_cases"] += 1
                                if first_remote_point != second_remote_point:
                                    examples.append((
                                        mode, n, seed, shared,
                                        "endpoint_closure",
                                        first_face, second_face,
                                        first_endpoint, second_endpoint,
                                    ))

                        for matching in (
                                ((first_pair[0], second_pair[0]),
                                 (first_pair[1], second_pair[1])),
                                ((first_pair[0], second_pair[1]),
                                 (first_pair[1], second_pair[0]))):
                            first_points = [
                                intersection(lines[shared], lines[left])
                                for left, _ in matching
                            ]
                            second_points = [
                                intersection(lines[shared], lines[right])
                                for _, right in matching
                            ]
                            if first_points != second_points:
                                continue
                            counts["shared_segment_face_pairs"] += 1
                            first_apex = intersection(
                                lines[first_pair[0]], lines[first_pair[1]])
                            second_apex = intersection(
                                lines[second_pair[0]], lines[second_pair[1]])
                            slope, intercept = lines[shared]
                            first_side = (
                                first_apex[1]
                                - slope * first_apex[0] - intercept)
                            second_side = (
                                second_apex[1]
                                - slope * second_apex[0] - intercept)
                            if first_side * second_side >= 0:
                                examples.append((
                                    mode, n, seed, shared,
                                    "shared_segment_side",
                                    first_face, second_face,
                                ))

    counts["violations"] = len(examples)
    result = {
        "status": "PASS" if not examples else "FAIL",
        "samples_per_mode": samples_per_mode,
        "n_range": [5, 10],
        "modes": [mode for mode, _, _ in MODES],
        "counts": counts,
        "sharpness_control": sharpness_control(),
        "violation_examples": examples[:20],
    }
    if examples:
        raise AssertionError(json.dumps(result, indent=2))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-per-mode", type=int, default=100)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = audit(args.samples_per_mode)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
