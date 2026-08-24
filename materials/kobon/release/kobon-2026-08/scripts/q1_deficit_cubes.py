#!/usr/bin/env python3
"""Build the exhaustive one-deficit cover of a canonical Q=1 gap-face target.

Requires sum(line capacities) - 3*T == 1.  Exact selection and the existing
per-line at-most constraints then imply that every model belongs to exactly
one cube: the unique line whose selected-side count is capacity minus one.
Each cube adds only that line's stricter at-most constraint.

Usage: q1_deficit_cubes.py N T OUT_DIR
"""
from __future__ import annotations

from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "math" / "kobon"))
sys.path.insert(0, str(ROOT / "scratch" / "kobon"))
import engine  # noqa: E402
import q1_gap_faces  # noqa: E402
from pysat.card import CardEnc, EncType  # noqa: E402


def build_cube(n, target, deficit_line):
    cnf, pool = q1_gap_faces.build(n, target)
    cap = n - 3 if deficit_line in (0, 1) else n - 2
    selected = [
        pool.id(("S",) + triple)
        for triple in combinations(range(n), 3)
        if deficit_line in triple
    ]
    cnf.extend(CardEnc.atmost(
        lits=selected, bound=cap - 1, vpool=pool,
        encoding=EncType.seqcounter,
    ))
    return cnf


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    n, target = int(sys.argv[1]), int(sys.argv[2])
    out_dir = Path(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)

    capacities = [n - 3, n - 3] + [n - 2] * (n - 2)
    slack = sum(capacities) - 3 * target
    if slack != 1:
        raise ValueError(f"one-deficit cover requires slack 1, got {slack}")

    rows = []
    for deficit_line in range(n):
        cnf = build_cube(n, target, deficit_line)
        path = out_dir / f"q1-deficit-{deficit_line:02d}.cnf"
        cnf.to_file(str(path))
        digest = sha256(path.read_bytes()).hexdigest()
        row = {
            "deficit_line": deficit_line,
            "capacity": capacities[deficit_line],
            "selected_sides": capacities[deficit_line] - 1,
            "vars": cnf.nv,
            "clauses": len(cnf.clauses),
            "bytes": path.stat().st_size,
            "sha256": digest,
            "cnf": str(path),
        }
        rows.append(row)
        print(json.dumps(row), flush=True)

    index = {
        "formulation": "canonical-q1-gap-faces-one-deficit-cover",
        "n": n,
        "target": target,
        "capacity_sum": sum(capacities),
        "selected_side_sum": 3 * target,
        "cover_argument": (
            "per-line counts are at most capacities; exact T counts each "
            "selected triple on exactly three lines; total slack is one"
        ),
        "cubes": rows,
    }
    index_path = out_dir / "index.json"
    index_path.write_text(json.dumps(index, indent=2) + "\n")
    print(json.dumps({"stage": "written", "index": str(index_path),
                      "cubes": len(rows)}), flush=True)


if __name__ == "__main__":
    main()
