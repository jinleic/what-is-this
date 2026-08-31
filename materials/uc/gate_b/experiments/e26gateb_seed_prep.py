#!/usr/bin/env python3
"""Full-set-pinned defect seeding for the e26 Gate B low-defect sweep.

The ``--score min-a`` walk in ``search_subtwofifths_q.py`` requires every seed
family to contain the full set ``[n]``, but the closure-defect minimizer in
``search_local_defect.py`` has no such constraint, so its frontier witnesses
usually cannot seed the walk.  This driver minimizes the same exact defect
count with one row pinned to the full set and never touched.  It writes
append-only, restart-indexed, deterministic schema-3 records in exactly the
pool format ``search_subtwofifths_q.seed_pool`` consumes; extra ``generator``
fields record the true origin, and the consumer re-verifies distinct rows, the
full set, cap, Reimer, and the exact defect count before use.  No float value
is persisted for promotion; checkpoints contain no timing fields.

Run from math/:
    OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
      uc/gate_b/experiments/e26gateb_seed_prep.py \
      --dimension 10 --sizes 13,15 --restarts 20 --steps 12000 --keep 24 \
      --seed 20260830 \
      --source-checkpoint uc/gate_b/experiments/local_defect_n10_checkpoint.jsonl \
      --checkpoint uc/gate_b/experiments/e26gateb_n10_fullsetseeds_checkpoint.jsonl
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

from search_local_defect import (  # noqa: E402
    coordinate_counts,
    defect_count,
    pool_families,
    read_witnesses,
    reimer_threshold,
    normalized,
    restart_seed,
)

SCHEMA = 3
GENERATOR = "uc/gate_b/experiments/e26gateb_seed_prep.py"
GENERATOR_NOTE = (
    "full-set-pinned closure-defect annealing for the e26 min-a sweep; "
    "search_subtwofifths_q.seed_pool re-verifies distinct rows, the full "
    "set, cap, Reimer, and the exact defect count before use"
)


def pinned_admissible(rows: list[int], dimension: int, size: int) -> bool:
    counts = coordinate_counts(tuple(rows), dimension)
    return (
        max(counts) <= 2 * size // 5
        and sum(counts) >= reimer_threshold(size)
    )

def greedy_repair(
    rows: list[int], dimension: int, rng: Random, max_tries: int = 40,
    attempts: int = 12,
) -> list[int] | None:
    """Walk rows containing [n] to admissibility by union-closure repair.

    Repeatedly replace a random row by a missing pairwise union until the
    family is union-closed, then fix the column cap by dropping the most
    loaded column's bit (re-closing after each drop) and raise incidence into
    a sparse column when Reimer needs slack.  Each pass re-closes, so the
    full-set row survives; returns None if a pass dead-ends.
    """
    rows_original = list(rows)
    size = len(rows)
    cap = 2 * size // 5
    need = reimer_threshold(size)
    full = (1 << dimension) - 1

    for _attempt in range(attempts):
        rows = list(rows_original)

        def close() -> None:
            for _ in range(max_tries * size):
                members = set(rows)
                holes = [
                    left | right
                    for left in rows for right in rows
                    if (left | right) not in members
                ]
                if not holes:
                    return
                union = rng.choice(holes)
                movable = [
                    index for index, row in enumerate(rows) if row != full
                ]
                rows[rng.choice(movable or [rng.randrange(size)])] = union

        close()
        for _pass in range(16 * max_tries):
            counts = coordinate_counts(tuple(rows), dimension)
            if max(counts) <= cap and sum(counts) >= need:
                if full not in set(rows):
                    return None
                return rows
            candidates = []
            for coordinate in range(dimension):
                if counts[coordinate] > cap:
                    holders = [
                        index for index, row in enumerate(rows)
                        if row >> coordinate & 1
                    ]
                    for index in holders:
                        proposed = rows[index] & ~(1 << coordinate)
                        if proposed not in set(rows):
                            candidates.append((index, proposed))
            if candidates:
                index, proposed = rng.choice(candidates)
                rows[index] = proposed
            else:
                spare = [
                    coordinate
                    for coordinate in range(dimension)
                    if counts[coordinate] < cap
                ]
                mixed = []
                for coordinate in spare:
                    missing = [
                        index for index, row in enumerate(rows)
                        if not row >> coordinate & 1
                    ]
                    for index in missing:
                        proposed = rows[index] | (1 << coordinate)
                        if proposed not in set(rows):
                            mixed.append((index, proposed))
                if not mixed:
                    break
                index, proposed = rng.choice(mixed)
                rows[index] = proposed
            close()
    return None


def full_set_p0(
    source_rows: list[int], size: int, dimension: int, rng: Random
) -> tuple[int, ...] | None:
    """A family of the given size containing [n]: derived from a source
    family by unioning rows, or None if the source cannot supply one."""
    source_size = len(source_rows)
    if source_size < size:
        return None
    rows = list(source_rows[:size])
    full = (1 << dimension) - 1
    if full not in set(rows):
        # Replace a row by the full set; if the replaced value already occurs
        # elsewhere, first swap that slot for a fresh unused row value.
        replace = rng.randrange(size)
        if rows[replace] in set(rows[:replace] + rows[replace + 1:]):
            pool = sorted(set(range(1 << dimension)) - set(rows) - {full})
            if not pool:
                return None
            rows[replace] = rng.choice(pool)
        rows[replace] = full
        if len(set(rows)) != size:
            return None
    fixed = greedy_repair(rows, dimension, rng)
    if fixed is None:
        return None
    return tuple(fixed)


def random_pinned_start(
    size: int, dimension: int, rng: Random
) -> tuple[int, ...] | None:
    """Distinct random rows with the full set pinned; walked to admissibility."""
    universe = 1 << dimension
    cap = 2 * size // 5
    need = reimer_threshold(size)
    if size > universe or need > dimension * cap:
        return None
    for _attempt in range(200):
        rows = rng.sample(range(universe - 1), size - 1)
        rows.append(universe - 1)
        movable = [index for index in range(size) if rows[index] != universe - 1]
        for _fix in range(64 * size):
            counts = coordinate_counts(tuple(rows), dimension)
            if max(counts) <= cap and sum(counts) >= need:
                return tuple(rows)
            if max(counts) > cap:
                worst = max(range(dimension), key=lambda c: counts[c])
                holders = [
                    index for index in movable if rows[index] >> worst & 1
                ]
                if not holders:
                    break
                index = rng.choice(holders)
                proposed = rows[index] & ~(1 << worst)
            else:
                spare = [c for c in range(dimension) if counts[c] < cap]
                if not spare:
                    break
                coordinate = rng.choice(spare)
                missing = [
                    index
                    for index in movable
                    if not rows[index] >> coordinate & 1
                ]
                if not missing:
                    break
                index = rng.choice(missing)
                proposed = rows[index] | (1 << coordinate)
            if proposed in rows:
                continue
            rows[index] = proposed
    return None


def anneal_pinned(
    size: int,
    dimension: int,
    degrees: tuple[int, ...],
    rng: Random,
    steps: int,
    keep: int,
    start_rows: tuple[int, ...] | None = None,
) -> dict[str, object]:
    """Minimize the exact defect count with the full-set row pinned."""
    full = (1 << dimension) - 1
    initial = None
    if (
        start_rows is not None
        and len(set(start_rows)) == size
        and full in set(start_rows)
    ):
        if pinned_admissible(list(start_rows), dimension, size):
            initial = tuple(start_rows)
    if initial is None:
        initial = random_pinned_start(size, dimension, rng)
    rows = list(initial or ())
    if not rows:
        return {"seeded": False}
    fixed = rows.index(full)
    movable = [index for index in range(size) if index != fixed]
    current = list(rows)
    current_bad = defect_count(tuple(current))
    best_bad = current_bad
    best_rows = tuple(sorted(current))
    pool: dict[tuple[int, ...], int] = {best_rows: best_bad}
    accepted = 0
    evaluated = 1
    hot = max(1.0, size ** 2 / 32.0)
    source_cap = 2 * size // 5
    source_need = reimer_threshold(size)

    def admissible_rows(candidate: list[int]) -> bool:
        if len(set(candidate)) != size:
            return False
        counts = coordinate_counts(tuple(candidate), dimension)
        return (
            max(counts) <= source_cap
            and sum(counts) >= source_need
        )

    repairs = 0
    for step in range(steps):
        temperature = hot * (1.0 - step / steps) + 1e-9
        draw = rng.random()
        previous = list(current)
        if draw < 0.20:
            victim = rng.choice(movable)
            current[victim] = rng.randrange(full)
        else:
            coordinate = rng.randrange(dimension)
            bit = 1 << coordinate
            ones = [
                index for index in movable
                if current[index] & bit and current[index] != full
            ]
            zeros = [
                index for index in movable if not current[index] & bit
            ]
            if not ones or not zeros:
                continue
            source = rng.choice(ones)
            target = rng.choice(zeros)
            current[source] ^= bit
            current[target] |= bit
        if not admissible_rows(current):
            current[:] = previous
            continue
        candidate_bad = defect_count(tuple(current))
        evaluated += 1
        delta = candidate_bad - current_bad
        if delta <= 0 or rng.random() < 2.0 ** (-delta / temperature):
            current_bad = candidate_bad
            accepted += 1
            if candidate_bad <= best_bad + max(2, size // 4):
                key = tuple(sorted(current))
                if key not in pool:
                    pool[key] = candidate_bad
            if candidate_bad < best_bad:
                best_bad = candidate_bad
                best_rows = tuple(sorted(current))

    ranked = sorted(pool.items(), key=lambda item: (item[1], item[0]))[:keep]
    return {
        "seeded": True,
        "degrees": list(degrees),
        "steps": steps,
        "evaluated": evaluated,
        "accepted": accepted,
        "best_defect_count": best_bad,
        "best_defect": str(Fraction(best_bad, size ** 2)),
        "best_family": list(best_rows),
        "best_normalized": normalized(best_rows, dimension),
        "pool": [
            {"family": list(family), "defect_count": bad}
            for family, bad in ranked
        ],
    }


def read_done(
    checkpoint: Path,
) -> dict[tuple[int, int, int], dict[str, object]]:
    done: dict[tuple[int, int, int], dict[str, object]] = {}
    if not checkpoint.exists():
        return done
    with checkpoint.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            if payload.get("schema") != SCHEMA:
                continue
            done[(payload["dimension"], payload["size"], payload["restart"])] = payload
    return done


def append_record(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--sizes", type=str, required=True)
    parser.add_argument("--restarts", type=int, default=16)
    parser.add_argument("--steps", type=int, default=12000)
    parser.add_argument("--keep", type=int, default=24)
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument("--source-checkpoint", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    args = parser.parse_args()

    sizes = tuple(int(token) for token in args.sizes.split(","))
    done = read_done(args.checkpoint)
    best_by_size: dict[str, object] = {}
    for size in sizes:
        restarts_done = 0
        restarted_from = 0
        best_bad = None
        best_family = None
        for restart in range(args.restarts):
            key = (args.dimension, size, restart)
            if key in done:
                restarted_from += 1
                payload = done[key]
            else:
                rng = Random(restart_seed(args.seed, args.dimension, size, restart))
                outcome = anneal_pinned(
                    size, args.dimension, (), rng, args.steps, args.keep
                )
                if not outcome.pop("seeded"):
                    outcome = {"best_defect_count": None}
                payload = {
                    "schema": SCHEMA,
                    "dimension": args.dimension,
                    "size": size,
                    "restart": restart,
                    "seed": args.seed,
                    "generator": GENERATOR,
                    "generator_note": GENERATOR_NOTE,
                    **outcome,
                }
                append_record(args.checkpoint, payload)
                done[key] = payload
            bad = payload.get("best_defect_count")
            if bad is not None and (best_bad is None or bad < best_bad):
                best_bad = bad
                best_family = payload.get("best_family")
            restarts_done += 1
        best_by_size[str(size)] = {
            "dimension": args.dimension,
            "size": size,
            "restarts_done": restarts_done,
            "restarted_from_checkpoint": restarted_from,
            "best_defect_count": best_bad,
            "best_defect": (
                str(Fraction(best_bad, size ** 2)) if best_bad is not None else None
            ),
            "best_family": best_family,
        }
    report = {
        "mode": "e26_gateb_fullset_seed_prep",
        "dimension": args.dimension,
        "sizes": list(sizes),
        "restarts": args.restarts,
        "steps": args.steps,
        "seed": args.seed,
        "generator": GENERATOR,
        "checkpoint": str(args.checkpoint),
        "best_by_size": best_by_size,
    }
    encoded = json.dumps(report, indent=2, sort_keys=True)
    sys.stdout.write(encoded + "\n")
    report_path = args.checkpoint.with_suffix(".seed_report.json")
    report_path.write_text(encoded + "\n")


if __name__ == "__main__":
    main()
