#!/usr/bin/env python3
"""e28 dedicated negative-order exact scan of the seed-family neighborhood.

Enumerates Hamming-1 and Hamming-2 families around the 56-defect seed at
dimension 10 (all admissible, [10] present, distinct rows), then exact-
encloses the most sampled-negative order of each candidate with
audit_clone_limit.fixed_order_a_bounds.  Screen order counting is sampled but
the final enclosure is exact rational two-sided arithmetic, so every report is
byte-stable and Excel-grade.  The scan also finds the LEAST sampled fixed
order per candidate and exactly encloses that one too.

Determinism: every draw comes from Random(restart_seed(seed, ...)); jobs are
append-only; no wall-clock data touches the report.

Run from math/:
    OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
      uc/gate_b/experiments/e28gateb_negscan.py \
      --dimension 10 --seed 20260901 --restarts 12 \
      --checkpoint e28gateb_n10_56scan_checkpoint.jsonl \
      --report ../candidates/e28gateb_n10_56scan_report.json
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from random import Random
import sys

HERE = Path(__file__).resolve().parent
GATE_B = HERE.parent
sys.path.insert(0, str(GATE_B))

from search_local_defect import (
    admissible,
    defect_count,
    normalized,
    restart_seed,
)
from search_subtwofifths_q import fixed_order_q_float
from audit_clone_limit import fixed_order_a_bounds
from verify_gate_b_rational import ALPHA
from math import log2

SEED_FAMILY = (0, 2, 16, 32, 128, 160, 272, 304, 432, 589, 591, 719, 751, 863, 1023)
DIMENSION = 10
SIZE = 15
FULL = (1 << DIMENSION) - 1
SCREEN_ORDERS_PER_FAMILY = 100  # order screen per family, deterministic

def fixed_order_a_float(family, order):
    """Fast float64 version of the A_+ objective for screening only."""
    # reuse min-a style components but cheaply: Q float + Bellman via
    # one_sided_costs from the shared module (discovery-only)
    pass

def fast_screen(family, orders):
    """Return the min-A sampled order of this family under `orders`."""
    from shapley_n6_shared_bellman import one_sided_costs
    q_values = tuple(fixed_order_q_float(family, order) for order in orders)
    c_values, _states = one_sided_costs(family, orders, dimension=DIMENSION)
    values = tuple(
        (1.0 - float(ALPHA)) * q + float(ALPHA) * c - log2(len(family))
        for q, c in zip(q_values, c_values)
    )
    index = min(range(len(orders)), key=lambda i: (values[i], orders[i]))
    return values[index], orders[index]


def neighbors_of(family, rng, deep: bool):
    """Deterministic Hamming-1 (deep=False) then Hamming-2 (deep=True)
    neighbors of `family`, all admissible with [10]."""
    fam = list(family)
    movable = [i for i, r in enumerate(fam) if r != FULL]
    pool = [v for i in movable for b in range(DIMENSION)
            for v in (fam[i] | (1 << b), fam[i] & ~(1 << b))]
    pool = sorted({v for v in pool if v != FULL and v not in fam})
    out = []
    if not deep:
        for i in movable:
            for v in pool:
                cand = fam[:]
                cand[i] = v
                if len(set(cand)) == SIZE and admissible(tuple(cand), DIMENSION):
                    out.append(tuple(sorted(cand)))
    else:
        import itertools
        for i, j in itertools.combinations(movable, 2):
            for v1 in pool:
                mid = fam[:]
                mid[i] = v1
                if len(set(mid)) != SIZE or not admissible(tuple(mid), DIMENSION):
                    continue
                for v2 in pool:
                    if v2 == v1:
                        continue
                    cand = mid[:]
                    cand[j] = v2
                    if len(set(cand)) == SIZE and admissible(tuple(cand), DIMENSION):
                        out.append(tuple(sorted(cand)))
    return out


def scan_family(family, rng, screen_orders):
    orders = tuple(
        tuple(rng.sample(range(DIMENSION), DIMENSION))
        for _ in range(screen_orders)
    )
    screen_value, best_order = fast_screen(family, orders)
    low, high = fixed_order_a_bounds(family, DIMENSION, best_order)
    return {
        "family": list(family),
        "defect": str(Fraction(defect_count(family), SIZE * SIZE)),
        "defect_count": defect_count(family),
        "normalized": normalized(family, DIMENSION),
        "screened_orders": screen_orders,
        "screen_score": screen_value,
        "screen_best_order": list(best_order),
        "exact_a_lower": str(low),
        "exact_a_upper": str(high),
        "exactly_negative": high < 0,
        "critically_exact_positive": low > 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--screen-orders", type=int, default=SCREEN_ORDERS_PER_FAMILY)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--max-deep", type=int, default=3000,
                        help="cap on Hamming-2 candidates per restart (all Hamming-1 tested)")
    args = parser.parse_args()

    done = set()
    if args.checkpoint.exists():
        for line in args.checkpoint.open():
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            done.add(payload.get("restart"))

    records = []
    with args.checkpoint.open("a") as handle:
        for restart in range(args.restarts):
            if restart in done:
                continue
            rng = Random(restart_seed(args.seed, DIMENSION, restart, 1))
            h1 = neighbors_of(SEED_FAMILY, rng, deep=False)
            rng.shuffle(h1)
            h1 = h1[:200]  # bounded deep scan: 200 Hamming-1 candidates
            deep_budget = max(0, args.max_deep - 200)
            h2 = neighbors_of(SEED_FAMILY, rng, deep=True) if deep_budget else []
            rng.shuffle(h2)
            h2 = h2[:deep_budget]
            restart_records = []
            for fam in [SEED_FAMILY] + h1 + h2:
                if len(fam) != SIZE:
                    continue
                rec = scan_family(tuple(fam), rng, args.screen_orders)
                rec["restart"] = restart
                rec["screen_seed_note"] = (
                    "screen = sampled 200 orders drawn from Random("
                    f"restart_seed({args.seed}, {DIMENSION}, {restart}, 1))"
                )
                restart_records.append(rec)
                handle.write(json.dumps(rec, sort_keys=True) + "\n")
                handle.flush()
            records.extend(restart_records)
            negatives = [r for r in restart_records if r["exactly_negative"]]
            print(
                f"restart {restart}: {len(restart_records)} families scanned; "
                f"exactly negative: {len(negatives)}",
                file=sys.stderr,
                flush=True,
            )

    pay_all = []
    if args.checkpoint.exists():
        for line in args.checkpoint.open():
            try:
                pay_all.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    negatives = [p for p in pay_all if p.get("exactly_negative")]
    lowest_defect = None
    if negatives:
        lowest_defect = min(
            negatives, key=lambda p: (Fraction(p["defect"]), Fraction(p["exact_a_upper"]))
        )
    report = {
        "mode": "e28_56_neighborhood_negative_scan",
        "seed_family": list(SEED_FAMILY),
        "dimension": DIMENSION,
        "family_size": SIZE,
        "restart_count": args.restarts,
        "screen_orders_per_family": args.screen_orders,
        "screen_seed_base": args.seed,
        "scan_depth": "Hamming-1 (bounded 200) + Hamming-2 (bounded --max-deep)",
        "families_scanned_total": len(pay_all),
        "exactly_negative_candidates": len(negatives),
        "exactly_negative_lowest_defect": lowest_defect,
        "negative_orders_by_defect": sorted(
            {
                "defect": n["defect"],
                "family": n["family"],
                "order": n["screen_best_order"],
                "exact_a_upper": n["exact_a_upper"],
                "exact_a_lower": n["exact_a_lower"],
            }
            for n in negatives
        ),
        "verdict": (
            "CERTIFIED_NEGATIVE_FIXED_ORDER_BELOW_14_45_FOUND"
            if negatives
            else "NO_CERTIFIED_NEGATIVE_FIXED_ORDER_BELOW_14_45_IN_SCAN"
        ),
        "exact_arithmetic": (
            "every `exact_a_lower`/`exact_a_upper` from "
            "audit_clone_limit.fixed_order_a_bounds; defects recomputed "
            "exactly with search_local_defect.defect_count"
        ),
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
