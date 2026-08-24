#!/usr/bin/env python3
"""Exact-rational audit of the selected-face vertex-sector bounds.

The audit samples deterministic arrangements with generic, parallel, and
concurrent degeneracies.  It checks the three local facts encoded by
engine.add_selected_face_sector_bounds:

* at most four triangular faces use any fixed pair of lines;
* at a finite multipoint, at most two use a fixed pair;
* a pair separated on both arcs of the cyclic slope order supports no face.

Finite sampling checks the implementation and label-order convention; the
universal justification is the local sector count at a line intersection.
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


def audit(samples_per_mode: int) -> dict:
    counts = {
        "arrangements": 0,
        "triangular_faces": 0,
        "line_pair_checks": 0,
        "multipoint_pair_checks": 0,
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
                    codegree = sum(i in face and j in face for face in faces)
                    counts["line_pair_checks"] += 1
                    if codegree > 4:
                        examples.append((mode, n, seed, i, j, "simple", codegree))

                    point = intersection(lines[i], lines[j])
                    if point is None:
                        continue
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
                            examples.append((mode, n, seed, i, j, "nonadjacent", codegree))

    counts["violations"] = len(examples)
    result = {
        "status": "PASS" if not examples else "FAIL",
        "samples_per_mode": samples_per_mode,
        "n_range": [5, 10],
        "modes": [mode for mode, _, _ in MODES],
        "counts": counts,
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
