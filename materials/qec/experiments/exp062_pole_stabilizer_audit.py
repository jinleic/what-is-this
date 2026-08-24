"""EXP-062: pole-stabilizer audit for the N=105 frontier.

The physical pole generator's translation stabilizer H controls two exact
reductions: pole supports are unions of H-cosets, and detector-sector anchors
need one coordinate per H-orbit.  This experiment tests whether H predicts the
observed easy [[210,18,8]] classes versus the six hard residuals.

Run:
  python experiments/exp062_pole_stabilizer_audit.py run
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E56 = _load("exp056_for_exp062", "exp056_odd_distance.py")
SHARD = ROOT / "results" / "partial_runs" / "exp055_screen" / "15x7.json"
OUT = ROOT / "results" / "processed" / "exp062_n105_pole_stabilizers.json"
SCHEMA = "exp062-pole-stabilizer-audit-v1"


def run() -> int:
    shard = json.loads(SHARD.read_text(encoding="utf-8"))
    rows = []
    for group in ("survivors", "undecided", "no_reference"):
        for index, record in enumerate(shard[group]):
            problem = E56.build_problem({
                "ell": record["ell"], "m": record["m"],
                "expected_k": record["k_parent"],
                "A": record["A"], "B": record["B"],
            })
            cover = E56.build_orbit_cover(problem)
            stabilizers = {sector["stabilizer_size"] for sector in cover["sectors"]}
            if len(stabilizers) != 1:
                raise RuntimeError("sector pole stabilizers disagree")
            H = stabilizers.pop()
            row = {
                "group": group, "index": index,
                "A": record["A"], "B": record["B"],
                "k": problem["k"], "threshold": record["threshold"],
                "witness_bound": record.get("witness_bound"),
                "pole_dimension": cover["pole_dimension"],
                "stabilizer_size": H,
                "coordinate_orbits_per_block": cover["sectors"][0]["anchor_orbit_count_per_block"],
                "pole_weights_divisible_by_H": True,
                "cover_complete": cover["complete"],
            }
            rows.append(row)
            print(json.dumps({k: row[k] for k in (
                "group", "index", "k", "stabilizer_size",
                "coordinate_orbits_per_block", "witness_bound",
            )}), flush=True)
    histogram = {}
    for row in rows:
        key = f"{row['group']}:H={row['stabilizer_size']}"
        histogram[key] = histogram.get(key, 0) + 1
    payload = {
        "schema": SCHEMA,
        "rows": rows,
        "histogram": histogram,
        "all_covers_complete": all(row["cover_complete"] for row in rows),
        "prediction": {
            "easy_survivors_large_H": all(
                row["stabilizer_size"] > 1 for row in rows if row["group"] == "survivors"
            ),
            "hard_residuals_H1": all(
                row["stabilizer_size"] == 1 for row in rows if row["group"] == "undecided"
            ),
        },
    }
    E56.atomic_write_json(OUT, payload)
    print(json.dumps({"histogram": histogram, "prediction": payload["prediction"]}, indent=1))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run",))
    parser.parse_args()
    return run()


if __name__ == "__main__":
    raise SystemExit(main())
