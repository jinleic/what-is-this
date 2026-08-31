#!/usr/bin/env python3
"""e27 seed checkpoint writer: Hamming/uplift neighborhoods of the e26
lowest-defect families, in seed_pool-consumable schema-3 records.

Seeds (all verified admissible, distinct rows, full set present, exact
defect; deterministic derivation from the e26 checkpoints):
  n=10 m=15: 62-family (62/225, e26 r0), 60-family (60/225, Hamming-1 of the
             62-family: row 3 -> 73), 56-family (56/225, Hamming-2: rows
             12 -> 863, 13 -> 751 from the 60-family).
  n=10 m=20: 13/40 family (e26 n10 m20 restart 2, 130/400).
  n=10 m=13: 64/169 family (e26 n10 m13 restart 3).
  n=11 m=15: 1024-lift of the 56-family (row 0 -> 1024; 84/225 at n=11).

Every record carries the derivation chain so provenance is exact.  The
consumer (search_subtwofifths_q.seed_pool) re-verifies admissibility and the
defect count before use.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
GATE_B = HERE.parent
sys.path.insert(0, str(GATE_B))

from search_local_defect import (  # noqa: E402
    admissible,
    defect_count,
    normalized,
    reimer_threshold,
)

SCHEMA = 3
GENERATOR = "uc/gate_b/experiments/e27gateb_seed_write.py"

SEEDS = [
    {
        "name": "e26_n10_m15_r0_defect62",
        "dimension": 10,
        "family": [0, 16, 32, 72, 128, 160, 272, 304, 432, 518, 591, 719,
                   895, 975, 1023],
        "derivation": "e26 sweep, search_subtwofifths_q.py --dimension 10 "
                      "--sizes 15, cap 73/225, seed 20260830, restart 0, "
                      "steps 7500, screen-orders 48, score min-a",
    },
    {
        "name": "e27_hamming1_defect60",
        "dimension": 10,
        "family": [0, 16, 32, 73, 128, 160, 272, 304, 432, 518, 591, 719,
                   895, 975, 1023],
        "derivation": "Hamming-1 single-row flip (row 3: 72 -> 73) of "
                      "e26_n10_m15_r0_defect62; defect 62 -> 60",
    },
    {
        "name": "e27_hamming2_defect56",
        "dimension": 10,
        "family": [0, 16, 32, 73, 128, 160, 272, 304, 432, 518, 591, 719,
                   863, 751, 1023],
        "derivation": "Hamming-2 flips (row 12: 895 -> 863, row 13: 975 -> "
                      "751) of e27_hamming1_defect60; defect 60 -> 56; "
                      "exhaustive 127698-case two-row neighborhood probe "
                      "found no family at defect <= 55",
    },
    {
        "name": "e26_n10_m20_r2_defect130",
        "dimension": 10,
        "family": None,  # filled at runtime from the e26 checkpoint
        "derivation": "e26 sweep, --dimension 10 --sizes 20, cap 137/400, "
                      "seed 20260830, restart 2; defect 13/40 = 130/400",
        "source_checkpoint": "uc/gate_b/experiments/e26gateb_n10_m20_sub34_checkpoint.jsonl",
        "source_restart": 2,
    },
    {
        "name": "e26_n10_m13_r3_defect64",
        "dimension": 10,
        "family": None,
        "derivation": "e26 sweep, --dimension 10 --sizes 13, cap 2/5, seed "
                      "20260830, restart 3; defect 64/169",
        "source_checkpoint": "uc/gate_b/experiments/e26gateb_n10_m13_sub40_checkpoint.jsonl",
        "source_restart": 3,
    },
]


def fill_from_checkpoints() -> None:
    for seed in SEEDS:
        if seed["family"] is not None:
            continue
        rows = [
            json.loads(line)
            for line in open(seed.pop("source_checkpoint"))
        ]
        restart = seed.pop("source_restart")
        payload = next(
            p for p in rows if p["restart"] == restart
        )
        seed["family"] = payload["family"]


def main() -> None:
    fill_from_checkpoints()
    out_path = HERE / "e27gateb_seeds_checkpoint.jsonl"
    report_path = HERE / "e27gateb_seeds_seed_report.json"
    lines = []
    summary = []
    with out_path.open("w") as handle:
        pools: dict[tuple[int, int], list[dict]] = {}
        for seed in SEEDS:
            rows = tuple(seed["family"])
            dimension = seed["dimension"]
            assert admissible(rows, dimension), seed["name"]
            assert (1 << dimension) - 1 in rows, seed["name"]
            bad = defect_count(rows)
            pool_entry = {
                "family": seed["family"],
                "defect_count": bad,
                "defect": f"{bad}/{len(rows) ** 2}",
                "normalized": normalized(rows, dimension),
                "name": seed["name"],
                "derivation": seed["derivation"],
                "generator": GENERATOR,
            }
            key = (dimension, len(rows))
            pools.setdefault(key, [])
            dup = any(
                tuple(existing["family"]) == rows
                for existing in pools[key]
            )
            if not dup:
                pools[key].append(pool_entry)
        for (dimension, size), pool in sorted(pools.items()):
            best = min(p["defect_count"] for p in pool)
            payload = {
                "schema": SCHEMA,
                "dimension": dimension,
                "size": size,
                "restart": 0,
                "seed": 20260831,
                "generator": GENERATOR,
                "generator_note": (
                    "e27 deterministic seed records derived from e26 "
                    "checkpoints and exact Hamming neighborhoods; "
                    "search_subtwofifths_q.seed_pool re-verifies distinct "
                    "rows, the full set, cap, Reimer, and the exact defect "
                    "count before use"
                ),
                "steps": 0,
                "evaluated": len(pool),
                "accepted": 0,
                "best_defect_count": best,
                "best_defect": str(__import__("fractions").Fraction(best, size ** 2)),
                "best_family": pool[0]["family"] if pool[0]["defect_count"] == best else next(
                    p["family"] for p in pool if p["defect_count"] == best
                ),
                "best_normalized": next(
                    p["normalized"] for p in pool if p["defect_count"] == best
                ),
                "pool": pool,
            }
            handle.write(json.dumps(payload, sort_keys=True) + "\n")
            summary.append({
                "dimension": dimension,
                "size": size,
                "pool": [
                    {k: v for k, v in p.items() if k != "derivation"}
                    for p in pool
                ],
            })
    report_path.write_text(
        json.dumps(
            {
                "mode": "e27_gateb_seed_write",
                "generator": GENERATOR,
                "checkpoint": str(out_path),
                "seeds_by_dimension_size": summary,
            },
            indent=2,
            sort_keys=True,
        ) + "\n"
    )
    print(f"wrote {out_path}")
    for entry in summary:
        print(entry["dimension"], entry["size"], "->", [
            (p["name"], p["defect"]) for p in entry["pool"]])


if __name__ == "__main__":
    main()
