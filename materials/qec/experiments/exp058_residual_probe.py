"""EXP-058 residual probe: decide the seven undecided [[186,10,d]] classes.

The EXP-055 screen ran its ladder reduced-witness -> CDCL (side z, 1e6
conflicts) -> exact CP-SAT (120 s) and left seven classes undecided.  The CP-SAT
route is a known dead end (EXP-056/next_actions).  This probe climbs BOTH sides
with a larger CDCL budget per cap: side x (never attempted before) first, then
side z with an escalating budget.  A SAT witness below the threshold dominates;
both-sides UNSAT at cap 13 proves survival and yields the exact distance.

Artifacts:
  results/partial_runs/exp058_residual_31x3.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT / "src"))

import exp055_odd_lattice_sweep as E55  # noqa: E402

from qec_research.distance.sat_decide import (  # noqa: E402
    css_side_instance,
    decide_weight_bounded,
)

SHARD = ROOT / "results" / "partial_runs" / "exp055_screen" / "31x3.json"
OUT = ROOT / "results" / "partial_runs" / "exp058_residual_31x3.json"
THRESHOLD = 12
SOLVER = "cadical195"


def probe_side(instance, cap: int, budget: int) -> dict[str, Any]:
    rec = decide_weight_bounded(
        instance, cap, solver_name=SOLVER, conflict_budget=budget
    )
    out = {
        "status": rec["status"],
        "cap": cap,
        "cnf_sha256": rec["cnf_sha256"],
        "solver": rec["solver"]["name"],
        "conflict_budget": budget,
        "wall_time_s": rec["solver"]["wall_time_s"],
        "stats": rec["solver"]["stats"],
    }
    if rec["status"] == "SAT":
        out["weight"] = rec["weight"]
        out["witness_support"] = [
            int(i) for i in np.flatnonzero(np.asarray(rec["vector"], dtype=np.uint8))
        ]
        out["verification"] = rec["verification"]
    return out


def main() -> int:
    shard = json.loads(SHARD.read_text())
    rows = shard["undecided"]
    if OUT.exists():
        done = json.loads(OUT.read_text())
    else:
        done = {"schema": "exp058-residual-31x3-v1", "records": []}
    finished = {json.dumps(r["B"]) for r in done["records"]}
    for row in rows:
        key = json.dumps(row["B"])
        if key in finished:
            continue
        ell, m = 31, 3
        HX, HZ = E55.E53.bb_from_terms(ell, m, row["A"], row["B"])
        instances = {
            side: css_side_instance(HX, HZ, side, block_length=ell * m)
            for side in ("x", "z")
        }
        import os
        only_side = os.environ.get("EXP058_ONLY_SIDE", "")
        solver_name = os.environ.get("EXP058_SOLVER", SOLVER)
        if only_side in ("x", "z"):
            instances = {only_side: instances[only_side]}
        record: dict[str, Any] = {
            "A": row["A"], "B": row["B"],
            "k_parent": row["k_parent"], "threshold": THRESHOLD,
            "reduced_witness_bound": row["witness_bound"],
            "probes": [],
        }
        # Side x first: never attempted by the screen.
        for side in ("x", "z"):
            if side not in instances:
                continue
            budget = int(os.environ.get(
                "EXP058_BUDGET",
                4_000_000 if side == "x" else 12_000_000,
            ))
            r12 = probe_side(instances[side], THRESHOLD, budget)
            r12["side"] = side
            record["probes"].append(r12)
            print(json.dumps({"B": row["B"], **r12}), flush=True)
            if r12["status"] == "SAT":
                record["verdict"] = "dominated"
                record["witness"] = {
                    "side": side, "weight": r12["weight"],
                    "support": r12["witness_support"],
                }
                break
        else:
            if all(p["status"] == "UNSAT" for p in record["probes"]):
                record["verdict"] = "lower_bound_13"
            else:
                record["verdict"] = "still_undecided"
        done["records"].append(record)
        OUT.write_text(json.dumps(done, indent=1, sort_keys=True) + "\n")
    print(json.dumps([r["verdict"] for r in done["records"]]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
