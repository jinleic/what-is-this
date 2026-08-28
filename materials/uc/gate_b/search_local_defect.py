#!/usr/bin/env python3
"""Minimum closure defect over cap-2/5 Reimer-admissible families.

Gate B refuted a universal scalar closure-defect repair with families whose
defect tends to one.  The surviving question is local: can an admissible family
be *nearly* union closed and still make ``A_+`` negative?  Its first half is
purely combinatorial and needs no entropy evaluation at all,

    eps_*(n, m) = min { eps_vee(F) : F admissible, |F| = m, F subseteq 2^[n] },

with ``eps_vee`` the ordered-pair-with-replacement closure defect of
``DEFINITIONS.md``.  This module computes that minimum.

Why the question is sharp
-------------------------
Admissibility forces every column degree into a narrow band.  At ``n = 6``,
``m = 25`` the cap is ``10`` and the exact integer Reimer threshold is ``59``,
while six columns at the cap supply only ``60``: the degree sequence is
``(10,10,10,10,10,10)`` or a permutation of ``(10,10,10,10,10,9)``.  Both Gate B
bases sit exactly on that boundary.  So the natural move is a *degree-preserving
swap*: move one set bit inside a column.  Every state it reaches is admissible by
construction, which is what makes the search efficient.

Two rigorous anchors
--------------------
1.  A family with ``eps_vee = 0`` is union closed, and admissibility caps every
    degree at ``floor(2m/5) < m/2``, so such a family would refute Frankl's
    conjecture.  Frankl is verified for ``n <= 12``, hence
    ``eps_*(n, m) > 0`` for every ``n <= 12`` computed here.
2.  Failing ordered pairs come in twos, because ``X | X = X`` never fails.  So a
    non-union-closed family has ``eps_vee >= 2/m^2``, and with ``-A_+ <= log2 m``
    the local ratio obeys ``-A_+/eps_vee <= m^2 log2(m) / 2``.  A local
    counterexample therefore needs ``m >= sqrt(2/eps)``: the local regime is
    reachable only through large families.

Modes
-----
``--mode exact`` enumerates *every* family of every admissible size on ``n``
coordinates and returns the exact minimum.  It is feasible for ``n <= 4``.

``--mode search`` runs restart-indexed simulated annealing on the defect count
alone.  Each restart is deterministic given ``(seed, n, m, degrees, index)``, and
one append-only JSONL record is written per completed restart, so a rerun skips
finished restarts and never rewrites history.  Emitted families are re-checked
for size, cap, exact integer incidence, distinct rows, and activity.

Coverage label: exact for ``n <= 4``; SEEDED HEURISTIC upper bounds on
``eps_*`` for larger ``n``.  Failure to find a smaller defect proves nothing.

Run with one thread and low priority:
    OMP_NUM_THREADS=1 nice -n 10 ./.venv/bin/python -B \
        uc/gate_b/search_local_defect.py --mode search --dimension 6
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import log2
from pathlib import Path
from random import Random
import sys
import time
from typing import Iterable, Iterator

HERE = Path(__file__).resolve().parent
DEFAULT_CHECKPOINT = HERE / "experiments" / "local_defect_checkpoint.jsonl"
DEFAULT_REPORT = HERE / "candidates" / "local_defect_frontier.json"
# 1: degree-preserving swaps, prescribed-degree seeding.
# 2: adds free row replacement and closure-repair proposals.
# 3: seeding walks to admissibility from distinct random rows, so large sizes
#    are reachable at all.  Earlier records are a different algorithm and are
#    counted as foreign rather than reused.
RECORD_SCHEMA = 3


# ---------------------------------------------------------------------------
# admissibility, all integer arithmetic
# ---------------------------------------------------------------------------


def reimer_threshold(size: int) -> int:
    """Least integer r with 2^(2r) >= size^size, i.e. ceil(m log2 m / 2)."""
    assert size >= 1
    target = size ** size
    low, high = 0, 1
    while (1 << (2 * high)) < target:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if (1 << (2 * middle)) >= target:
            high = middle
        else:
            low = middle + 1
    return low


def coordinate_counts(rows: Iterable[int], dimension: int) -> tuple[int, ...]:
    counts = [0] * dimension
    for row in rows:
        for coordinate in range(dimension):
            if (row >> coordinate) & 1:
                counts[coordinate] += 1
    return tuple(counts)


def admissible(rows: tuple[int, ...], dimension: int) -> bool:
    size = len(rows)
    if size < 3 or len(set(rows)) != size:
        return False
    if max(rows) >= (1 << dimension):
        return False
    counts = coordinate_counts(rows, dimension)
    return (
        max(counts) <= 2 * size // 5
        and sum(counts) >= reimer_threshold(size)
    )


def normalized(rows: tuple[int, ...], dimension: int) -> bool:
    columns = tuple(
        tuple((row >> coordinate) & 1 for row in rows)
        for coordinate in range(dimension)
    )
    return all(any(column) for column in columns) and len(set(columns)) == dimension


def feasible_sizes(dimension: int) -> tuple[int, ...]:
    return tuple(
        size for size in range(3, (1 << dimension) + 1)
        if reimer_threshold(size) <= dimension * (2 * size // 5)
    )


# ---------------------------------------------------------------------------
# closure defect
# ---------------------------------------------------------------------------


def defect_count(rows: tuple[int, ...]) -> int:
    """Number of ordered pairs (X, Y) in F^2 with X | Y not in F."""
    members = set(rows)
    bad = 0
    for left in rows:
        for right in rows:
            if (left | right) not in members:
                bad += 1
    return bad


def defect(rows: tuple[int, ...]) -> Fraction:
    return Fraction(defect_count(rows), len(rows) ** 2)


# ---------------------------------------------------------------------------
# exact enumeration
# ---------------------------------------------------------------------------


def exact_minimum(dimension: int) -> dict[str, object]:
    universe = 1 << dimension
    results = {}
    for size in feasible_sizes(dimension):
        cap = 2 * size // 5
        need = reimer_threshold(size)
        best_bad = None
        best_rows = None
        best_normalized_bad = None
        best_normalized_rows = None
        examined = 0
        for rows in combinations(range(universe), size):
            counts = coordinate_counts(rows, dimension)
            if max(counts) > cap or sum(counts) < need:
                continue
            examined += 1
            bad = defect_count(rows)
            if best_bad is None or bad < best_bad:
                best_bad, best_rows = bad, rows
            if normalized(rows, dimension):
                if best_normalized_bad is None or bad < best_normalized_bad:
                    best_normalized_bad, best_normalized_rows = bad, rows
        if best_bad is None:
            continue
        results[str(size)] = {
            "size": size,
            "cap": cap,
            "incidence_need": need,
            "admissible_count": examined,
            "min_defect_count": best_bad,
            "min_defect": str(Fraction(best_bad, size ** 2)),
            "min_defect_decimal": float(Fraction(best_bad, size ** 2)),
            "witness": list(best_rows),
            "min_defect_count_normalized": best_normalized_bad,
            "min_defect_normalized": (
                None if best_normalized_bad is None
                else str(Fraction(best_normalized_bad, size ** 2))
            ),
            "witness_normalized": (
                None if best_normalized_rows is None
                else list(best_normalized_rows)
            ),
        }
    return {
        "mode": "exact",
        "dimension": dimension,
        "sizes": results,
        "complete": True,
    }


# ---------------------------------------------------------------------------
# degree-preserving annealing


def restart_seed(seed: int, dimension: int, size: int, restart: int) -> int:
    """Deterministic integer seed; `Random` rejects tuples and `hash` is salted."""
    return (((seed * 1_000_003 + dimension) * 10_007 + size) * 101 + restart)


def degree_sequences(size: int, dimension: int, rng: Random) -> tuple[int, ...]:
    """Random feasible column-degree sequence, biased to the cap boundary."""
    cap = 2 * size // 5
    need = reimer_threshold(size)
    total = dimension * cap
    assert need <= total
    degrees = [cap] * dimension
    slack = total - need
    for _ in range(rng.randrange(slack + 1)):
        choices = [c for c in range(dimension) if degrees[c] > 0]
        if not choices:
            break
        degrees[rng.choice(choices)] -= 1
    assert sum(degrees) >= need and max(degrees) <= cap
    return tuple(degrees)


def seed_family(
    size: int, dimension: int, degrees: tuple[int, ...], rng: Random
) -> tuple[int, ...] | None:
    """Any admissible family of the given size; `degrees` is a soft target.

    Sampling column by column from a prescribed degree sequence collapses at
    large sizes: at n=7, m=45 the mean row weight is 124/45 < 3, so independent
    column draws almost always repeat a row and every attempt is rejected.  That
    silently dropped m=45 -- the size of the published base -- from the first
    sweep.  Instead start from distinct random rows and walk to admissibility by
    removing a set bit from the most loaded column, or adding one to a column
    under the cap, never creating a duplicate row.
    """
    universe = 1 << dimension
    if size > universe:
        return None
    cap = 2 * size // 5
    need = reimer_threshold(size)
    if need > dimension * cap:
        return None

    for _attempt in range(200):
        rows = rng.sample(range(universe), size)
        for _fix in range(64 * size):
            counts = coordinate_counts(rows, dimension)
            if max(counts) <= cap and sum(counts) >= need:
                return tuple(rows)
            present = set(rows)
            if max(counts) > cap:
                worst = max(range(dimension), key=lambda c: counts[c])
                holders = [i for i, row in enumerate(rows) if row >> worst & 1]
                index = rng.choice(holders)
                proposed = rows[index] & ~(1 << worst)
            else:
                spare = [c for c in range(dimension) if counts[c] < cap]
                if not spare:
                    break
                coordinate = rng.choice(spare)
                missing = [
                    i for i, row in enumerate(rows)
                    if not row >> coordinate & 1
                ]
                if not missing:
                    break
                index = rng.choice(missing)
                proposed = rows[index] | (1 << coordinate)
            if proposed in present:
                continue
            rows[index] = proposed
    return None


def anneal(
    size: int,
    dimension: int,
    degrees: tuple[int, ...],
    rng: Random,
    steps: int,
    keep: int,
) -> dict[str, object]:
    """Minimize the defect count.

    Three proposal kinds.  A degree-preserving swap moves one set bit inside a
    column, so it never leaves the admissible class and needs no re-check.  A
    free row replacement changes the degree sequence and is re-checked against
    cap and the exact integer incidence threshold.  A repair move picks an
    actually failing pair and tries to install its union as a row, which is the
    only proposal that directly attacks the objective rather than diffusing.
    """
    rows = seed_family(size, dimension, degrees, rng)
    if rows is None:
        return {"seeded": False}
    current = list(rows)
    current_bad = defect_count(tuple(current))
    best_bad = current_bad
    best_rows = tuple(sorted(current))
    pool: dict[tuple[int, ...], int] = {best_rows: best_bad}
    accepted = 0
    evaluated = 1
    hot = max(1.0, size ** 2 / 32.0)

    def admissible_rows(candidate: list[int]) -> bool:
        if len(set(candidate)) != len(candidate):
            return False
        counts = coordinate_counts(candidate, dimension)
        return (
            max(counts) <= 2 * size // 5
            and sum(counts) >= reimer_threshold(size)
        )

    def missing_union(candidate: list[int]) -> int | None:
        members = set(candidate)
        holes = [
            left | right
            for left in candidate for right in candidate
            if (left | right) not in members
        ]
        return rng.choice(holes) if holes else None

    repairs = 0
    for step in range(steps):
        temperature = hot * (1.0 - step / steps) + 1e-9
        draw = rng.random()
        previous = list(current)

        if draw < 0.15:
            union = missing_union(current)
            if union is None:
                break
            victim = rng.randrange(size)
            current[victim] = union
            if not admissible_rows(current):
                current[:] = previous
                continue
            repairs += 1
        elif draw < 0.30:
            victim = rng.randrange(size)
            current[victim] = rng.randrange(1 << dimension)
            if not admissible_rows(current):
                current[:] = previous
                continue
        else:
            coordinate = rng.randrange(dimension)
            bit = 1 << coordinate
            ones = [index for index, row in enumerate(current) if row & bit]
            zeros = [index for index, row in enumerate(current) if not row & bit]
            if not ones or not zeros:
                continue
            source = rng.choice(ones)
            target = rng.choice(zeros)
            current[source] ^= bit
            current[target] |= bit
            if len(set(current)) != size:
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
        else:
            current[:] = previous

    ranked = sorted(pool.items(), key=lambda item: (item[1], item[0]))[:keep]
    return {
        "seeded": True,
        "degrees": list(degrees),
        "steps": steps,
        "evaluated": evaluated,
        "accepted": accepted,
        "repairs": repairs,
        "best_defect_count": best_bad,
        "best_defect": str(Fraction(best_bad, size ** 2)),
        "best_family": list(best_rows),
        "best_normalized": normalized(best_rows, dimension),
        "pool": [
            {"family": list(family), "defect_count": bad}
            for family, bad in ranked
        ],
    }


def read_checkpoint(path: Path) -> tuple[dict[tuple[int, int, int], dict], int, int]:
    """Return (finished restarts, foreign-schema lines, truncated lines)."""
    done: dict[tuple[int, int, int], dict] = {}
    foreign = truncated = 0
    if not path.exists():
        return done, foreign, truncated
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                truncated += 1
                continue
            if payload.get("schema") != RECORD_SCHEMA:
                foreign += 1
                continue
            key = (payload["dimension"], payload["size"], payload["restart"])
            done[key] = payload
    return done, foreign, truncated


def append_record(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()


def read_witnesses(path: Path, dimension: int) -> dict[int, dict]:
    """Best exactly re-verified witness per size, across every record schema.

    Schema strictness governs *work skipping*: a restart finished by a different
    algorithm must not mark this algorithm's restart as done.  It must not
    govern *evidence*.  A witness is a concrete family, and every family here is
    re-checked for admissibility and its defect recount before it is allowed to
    set a bound, so an older algorithm's verified family is still a valid upper
    bound on eps_*.  The schema that produced it is recorded.
    """
    best: dict[int, dict] = {}
    if not path.exists():
        return best
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            if payload.get("dimension") != dimension:
                continue
            candidates = []
            family = payload.get("best_family")
            if family:
                candidates.append(tuple(family))
            for entry in payload.get("pool", []):
                candidates.append(tuple(entry["family"]))
            for rows in candidates:
                if not admissible(rows, dimension):
                    continue
                size = len(rows)
                bad = defect_count(rows)
                current = best.get(size)
                if current is None or bad < current["defect_count"]:
                    best[size] = {
                        "size": size,
                        "cap": 2 * size // 5,
                        "incidence_need": reimer_threshold(size),
                        "pair_floor": str(Fraction(2, size ** 2)),
                        "defect_count": bad,
                        "defect": str(Fraction(bad, size ** 2)),
                        "defect_decimal": float(Fraction(bad, size ** 2)),
                        "family": list(rows),
                        "normalized": normalized(rows, dimension),
                        "from_schema": payload.get("schema"),
                        "reverified": True,
                    }
    return best


def search(
    dimension: int,
    sizes: tuple[int, ...],
    restarts: int,
    steps: int,
    keep: int,
    seed: int,
    checkpoint: Path,
) -> dict[str, object]:
    done, foreign, truncated = read_checkpoint(checkpoint)
    started = time.monotonic()
    resumed = 0
    for size in sizes:
        for restart in range(restarts):
            key = (dimension, size, restart)
            if key in done:
                resumed += 1
                continue
            rng = Random(restart_seed(seed, dimension, size, restart))
            degrees = degree_sequences(size, dimension, rng)
            outcome = anneal(size, dimension, degrees, rng, steps, keep)
            if not outcome.pop("seeded"):
                outcome = {"degrees": list(degrees), "best_defect_count": None}
            family = outcome.get("best_family")
            if family is not None:
                rows = tuple(family)
                assert admissible(rows, dimension), (dimension, size, rows)
                assert defect_count(rows) == outcome["best_defect_count"]
            payload = {
                "schema": RECORD_SCHEMA,
                "dimension": dimension,
                "size": size,
                "restart": restart,
                "seed": seed,
                **outcome,
            }
            append_record(checkpoint, payload)
            done[key] = payload
            elapsed = time.monotonic() - started
            print(
                "n=%d m=%d restart %d/%d best_bad=%s (%.1f s)"
                % (dimension, size, restart + 1, restarts,
                   payload.get("best_defect_count"), elapsed),
                file=sys.stderr, flush=True,
            )

    frontier: dict[str, object] = {}
    for (dim, size, _restart), payload in sorted(done.items()):
        if dim != dimension or size not in sizes:
            continue
        bad = payload.get("best_defect_count")
        if bad is None:
            continue
        entry = frontier.setdefault(str(size), {
            "size": size,
            "cap": 2 * size // 5,
            "incidence_need": reimer_threshold(size),
            "pair_floor": str(Fraction(2, size ** 2)),
            "best_defect_count": bad,
            "best_defect": str(Fraction(bad, size ** 2)),
            "best_family": payload["best_family"],
            "best_normalized": payload.get("best_normalized"),
            "restarts": 0,
        })
        entry["restarts"] += 1
        if bad < entry["best_defect_count"]:
            entry.update(
                best_defect_count=bad,
                best_defect=str(Fraction(bad, size ** 2)),
                best_family=payload["best_family"],
                best_normalized=payload.get("best_normalized"),
            )
    return {
        "mode": "search",
        "dimension": dimension,
        "seed": seed,
        "restarts": restarts,
        "steps": steps,
        "resumed_restarts": resumed,
        "foreign_schema_lines": foreign,
        "truncated_checkpoint_lines": truncated,
        "runtime_seconds": round(time.monotonic() - started, 2),
        "sizes": frontier,
        "complete": False,
        "coverage": "SEEDED HEURISTIC upper bound on eps_*",
    }


def pool_families(checkpoint: Path, dimension: int, limit: int) -> list[dict]:
    """Lowest-defect distinct admissible families recorded so far."""
    done, _foreign, _truncated = read_checkpoint(checkpoint)
    seen: dict[tuple[int, ...], int] = {}
    for (dim, _size, _restart), payload in done.items():
        if dim != dimension:
            continue
        for entry in payload.get("pool", []):
            rows = tuple(entry["family"])
            bad = entry["defect_count"]
            if rows not in seen or bad < seen[rows]:
                seen[rows] = bad
    ranked = sorted(seen.items(), key=lambda item: (item[1] / len(item[0]) ** 2, item[1]))
    out = []
    for rows, bad in ranked[:limit]:
        assert admissible(rows, dimension)
        out.append({
            "family": list(rows),
            "size": len(rows),
            "defect_count": bad,
            "defect": str(Fraction(bad, len(rows) ** 2)),
            "normalized": normalized(rows, dimension),
        })
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("exact", "search", "pool", "frontier"),
        default="search",
    )
    parser.add_argument("--dimension", type=int, default=6)
    parser.add_argument("--sizes", type=str, default="")
    parser.add_argument("--restarts", type=int, default=6)
    parser.add_argument("--steps", type=int, default=4000)
    parser.add_argument("--keep", type=int, default=40)
    parser.add_argument("--limit", type=int, default=60)
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()

    if args.mode == "exact":
        result = exact_minimum(args.dimension)
    elif args.mode == "frontier":
        witnesses = read_witnesses(args.checkpoint, args.dimension)
        ranked = sorted(witnesses.values(), key=lambda row: row["defect_decimal"])
        result = {
            "mode": "frontier",
            "dimension": args.dimension,
            "sizes": {str(row["size"]): row for row in witnesses.values()},
            "min_defect_over_sizes": ranked[0] if ranked else None,
            "min_defect_at_size_at_least_20": min(
                (row for row in witnesses.values() if row["size"] >= 20),
                key=lambda row: row["defect_decimal"],
                default=None,
            ),
            "note": (
                "upper bounds on eps_*(n, m); every witness re-verified for "
                "admissibility and defect count, schema of origin recorded"
            ),
        }
    elif args.mode == "pool":
        result = {
            "mode": "pool",
            "dimension": args.dimension,
            "families": pool_families(args.checkpoint, args.dimension, args.limit),
        }
    else:
        if args.sizes:
            sizes = tuple(int(token) for token in args.sizes.split(","))
        else:
            sizes = feasible_sizes(args.dimension)
        result = search(
            args.dimension, sizes, args.restarts, args.steps,
            args.keep, args.seed, args.checkpoint,
        )

    encoded = json.dumps(result, indent=2, sort_keys=True)
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded + "\n")
    print(encoded)


if __name__ == "__main__":
    main()
