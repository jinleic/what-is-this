#!/usr/bin/env python3
"""Hunt admissible families with a large closure-defect repair ratio.

Gate B established

    c_cl^star = sup { -A_+(F) / eps_vee(F) } = +infinity

using Cartesian powers of two bases whose defect tends to *one*.  The surviving
question is local: can the ratio be large while the defect is *small*?  By
Lemma L2 (defect multiplies under products) the power mechanism can never
produce such a family, so the local regime needs a genuinely different object.

Search objective
----------------
For a coefficient ``lam`` the scalarization

    score(F) = A_+(F) + lam * eps_vee(F)

is exactly the proposed repair functional: ``score(F) < 0`` iff
``-A_+(F)/eps_vee(F) > lam``.  Minimizing ``score`` over admissible families is
therefore a direct attack on the coefficient ``lam``, and a sweep over ``lam``
traces the reachable ratio frontier at fixed dimension.  Because ``eps_vee``
enters with a positive weight, low-defect families are actively preferred: this
is the search formulation that the earlier ``A_+``-only censuses lacked.

Calibration.  The known finite bases give true ratios

    B (n=7, m=45): -A_+/eps_vee = 0.0286491.../ (64/81)  = 0.0362591270...
    D (n=6, m=25): -A_+/eps_vee = 0.0136721.../(444/625) = 0.0192456471...

so any family with a ratio above ``0.0362592`` is a *stronger finite base* than
the one in the published Gate B certificate, and any negative family with
defect below ``64/81`` moves the local frontier down.

Arithmetic and labels
---------------------
``eps_vee`` is exact (a ``Fraction`` over ordered pairs with replacement).
``A_+`` is float64 here, over a fixed sampled order set during the walk and over
*all* ``n!`` orders for the finalists.  Every reported family is re-checked for
size, cap, exact integer incidence, distinct rows, and activity.  Coverage is
SEEDED HEURISTIC: a negative find must still be certified by
``verify_gate_b_rational.py`` before it can carry a PROVED label, and failure to
find one proves nothing.

Resumability: one append-only JSONL record per ``(dimension, lam, size,
restart)`` with an explicit ``schema``.  A rerun skips finished restarts, and
records written by another schema are counted and ignored rather than trusted.

Run with one thread and low priority:
    OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
        uc/gate_b/hunt_local_ratio.py --dimension 7 --lam 0.05
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import permutations
import json
from math import log2
from pathlib import Path
from random import Random
import sys
import time

HERE = Path(__file__).resolve().parent
UC_ROOT = HERE.parent
sys.path.insert(0, str(UC_ROOT))
sys.path.insert(0, str(HERE))

from shapley_entropy import ALPHA  # noqa: E402
from shapley_join_loss import shapley_iid  # noqa: E402
from shapley_n6_shared_bellman import one_sided_costs  # noqa: E402

from search_local_defect import (  # noqa: E402
    admissible,
    coordinate_counts,
    defect_count,
    degree_sequences,
    normalized,
    reimer_threshold,
    restart_seed,
    seed_family,
)

DEFAULT_CHECKPOINT = HERE / "experiments" / "local_ratio_checkpoint.jsonl"
DEFAULT_REPORT = HERE / "candidates" / "local_ratio_frontier.json"
RECORD_SCHEMA = 1
SCREEN_ORDERS = 12


def objective(family: tuple[int, ...], dimension: int, orders) -> tuple[float, float, float]:
    """Return (A_+, Q, C_+) with C_+ averaged over the given orders."""
    iid = shapley_iid(family, dimension)
    costs, _states = one_sided_costs(family, orders, dimension=dimension)
    one_sided = sum(costs) / len(costs)
    value = (1.0 - ALPHA) * iid + ALPHA * one_sided - log2(len(family))
    return value, iid, one_sided


def full_objective(family: tuple[int, ...], dimension: int):
    orders = tuple(permutations(range(dimension)))
    return objective(family, dimension, orders)


def ratio(a_plus: float, defect: Fraction) -> float | None:
    if defect == 0:
        return None
    return -a_plus / float(defect)


def known_base(dimension: int) -> tuple[int, ...]:
    """The certified Gate B base for this dimension, as a feasible start."""
    if dimension == 6:
        cells = {(0, 0), (0, 1), (1, 0), (1, 2), (2, 1)}
        return tuple(sorted(
            left | (right << 3)
            for left in range(8) for right in range(8)
            if (bin(left).count("1"), bin(right).count("1")) in cells
        ))
    if dimension == 7:
        payload = json.loads(
            (HERE / "candidates" / "n7_block_extremizer.json").read_text())
        return tuple(payload["family_rows"])
    raise ValueError(f"no certified base for dimension {dimension}")


def seed_pool(
    paths: tuple[Path, ...], dimension: int, size: int
) -> list[tuple[tuple[int, ...], int, float | None]]:
    """Re-verified candidate seeds of the given size, with any known ``A_+``.

    A sweep seeded only from the minimum-defect witness cannot climb back into
    the large-defect negativity region: at cap 0.79 it reported ``+0.037`` while
    the certified base ``B`` (defect ``64/81``, ``A_+ = -0.0286``) was feasible.
    Seeding instead from the best family known under each cap makes every point
    of the curve at least as good as everything found so far.
    """
    seen: dict[tuple[int, ...], tuple[int, float | None]] = {}

    def offer(rows_list, known: float | None) -> None:
        rows = tuple(sorted(rows_list))
        if len(rows) != size or not admissible(rows, dimension):
            return
        bad = defect_count(rows)
        previous = seen.get(rows)
        if previous is None:
            seen[rows] = (bad, known)
        elif known is not None and (
            previous[1] is None or known < previous[1]
        ):
            seen[rows] = (bad, known)

    for path in paths:
        if not path.exists():
            continue
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
                for key in ("best_family", "family"):
                    if payload.get(key):
                        offer(payload[key], payload.get("A_plus"))
                for entry in payload.get("pool", []):
                    offer(entry["family"], None)
                for entry in payload.get("frontier", []):
                    offer(entry["family"], entry.get("A_plus"))
    return [(rows, bad, known) for rows, (bad, known) in seen.items()]


def best_seed(
    pool: list[tuple[tuple[int, ...], int, float | None]], budget: Fraction, size: int
) -> tuple[int, ...] | None:
    """Feasible seed under the defect cap: lowest known ``A_+``, else lowest defect."""
    feasible = [
        entry for entry in pool
        if Fraction(entry[1], size ** 2) <= budget
    ]
    if not feasible:
        return None
    scored = [entry for entry in feasible if entry[2] is not None]
    if scored:
        return min(scored, key=lambda entry: entry[2])[0]
    return min(feasible, key=lambda entry: entry[1])[0]


def sweep(
    dimension: int,
    size: int,
    defect_cap: Fraction,
    start: tuple[int, ...],
    rng: Random,
    steps: int,
    screen: int,
) -> dict[str, object]:
    """Minimize A_+ subject to a hard cap on the closure defect.

    This measures the local frontier function

        Lambda_m(e) = min { A_+(F) : F admissible, |F| = m, eps_vee(F) <= e },

    whose sign change locates
    ``e* = min { eps_vee(F) : F admissible, A_+(F) < 0 }`` -- the single number
    that decides whether the local regime is populated.  The defect cap is a
    hard filter rather than a penalty, so every visited family is feasible and
    the reported minimum is a genuine upper bound on ``Lambda_m(e)``.
    """
    cap = 2 * size // 5
    need = reimer_threshold(size)
    denominator = size ** 2
    budget = defect_cap * denominator
    orders = (
        tuple(permutations(range(dimension))) if dimension <= 6
        else tuple(
            tuple(rng.sample(range(dimension), dimension)) for _ in range(screen)
        )
    )

    def feasible(candidate: list[int]) -> bool:
        if len(set(candidate)) != len(candidate):
            return False
        counts = coordinate_counts(candidate, dimension)
        if max(counts) > cap or sum(counts) < need:
            return False
        return defect_count(tuple(candidate)) <= budget

    assert admissible(start, dimension) and defect_count(start) <= budget
    current = list(start)
    current_a = objective(tuple(sorted(current)), dimension, orders)[0]
    best_a, best_rows = current_a, tuple(sorted(current))
    evaluated = 1
    hot = 0.05

    for step in range(steps):
        temperature = hot * (1.0 - step / steps) + 1e-9
        previous = list(current)
        draw = rng.random()
        if draw < 0.20:
            current[rng.randrange(size)] = rng.randrange(1 << dimension)
        else:
            coordinate = rng.randrange(dimension)
            bit = 1 << coordinate
            ones = [index for index, row in enumerate(current) if row & bit]
            zeros = [index for index, row in enumerate(current) if not row & bit]
            if not ones or not zeros:
                continue
            current[rng.choice(ones)] ^= bit
            current[rng.choice(zeros)] |= bit
        if not feasible(current):
            current[:] = previous
            continue
        key = tuple(sorted(current))
        candidate_a = objective(key, dimension, orders)[0]
        evaluated += 1
        delta = candidate_a - current_a
        if delta <= 0 or rng.random() < 2.0 ** (-delta / temperature):
            current_a = candidate_a
            if candidate_a < best_a:
                best_a, best_rows = candidate_a, key
        else:
            current[:] = previous

    full_a, _q, _c = full_objective(best_rows, dimension)
    exact_defect = Fraction(defect_count(best_rows), denominator)
    assert admissible(best_rows, dimension)
    assert exact_defect <= defect_cap
    return {
        "size": size,
        "defect_cap": str(defect_cap),
        "defect_cap_decimal": float(defect_cap),
        "steps": steps,
        "evaluated": evaluated,
        "objective_orders": "all" if dimension <= 6 else screen,
        "min_A_plus_screened": best_a,
        "min_A_plus": full_a,
        "defect": str(exact_defect),
        "defect_decimal": float(exact_defect),
        "ratio": ratio(full_a, exact_defect),
        "negative": full_a < 0.0,
        "family": list(best_rows),
        "normalized": normalized(best_rows, dimension),
    }


def descend(
    dimension: int,
    rng: Random,
    steps: int,
    margin: float,
    penalty: float,
    screen: int,
) -> dict[str, object]:
    """Push the defect down from a certified negative family.

    The random-restart hunt above never reached a negative family at all, and in
    particular never rediscovered the two certified bases, so its failures
    measured search difficulty rather than the landscape.  This walk instead
    *starts* at a certified negative family and minimizes

        score(F) = eps_vee(F) + penalty * max(0, A_+(F) + margin),

    which leaves the objective equal to the defect exactly while
    ``A_+ <= -margin`` holds and charges heavily for losing negativity.  Every
    negative family it visits is a witness for the local frontier, so the walk
    reports the whole Pareto set of (defect, A_+) over negatives, not just its
    endpoint.
    """
    start = known_base(dimension)
    size = len(start)
    cap = 2 * size // 5
    need = reimer_threshold(size)
    denominator = size ** 2
    exact_orders = tuple(permutations(range(dimension)))
    use_full = dimension <= 6
    orders = exact_orders if use_full else tuple(
        tuple(rng.sample(range(dimension), dimension)) for _ in range(screen)
    )

    def admissible_rows(candidate: list[int]) -> bool:
        if len(set(candidate)) != len(candidate):
            return False
        counts = coordinate_counts(candidate, dimension)
        return max(counts) <= cap and sum(counts) >= need

    def score_of(rows: tuple[int, ...]) -> tuple[float, float, int]:
        bad = defect_count(rows)
        a_plus, _q, _c = objective(rows, dimension, orders)
        value = bad / denominator + penalty * max(0.0, a_plus + margin)
        return value, a_plus, bad

    assert admissible(start, dimension)
    current = list(start)
    current_score, current_a, current_bad = score_of(start)
    assert current_a < 0.0, "start must be a negative family"
    negatives: dict[tuple[int, ...], tuple[int, float]] = {
        start: (current_bad, current_a)
    }
    best = (current_bad, start, current_a)
    evaluated = 1
    accepted = 0
    hot = 0.02

    for step in range(steps):
        temperature = hot * (1.0 - step / steps) + 1e-9
        previous = list(current)
        draw = rng.random()
        if draw < 0.25:
            members = set(current)
            holes = [
                left | right
                for left in current for right in current
                if (left | right) not in members
            ]
            if not holes:
                break
            current[rng.randrange(size)] = rng.choice(holes)
        elif draw < 0.40:
            current[rng.randrange(size)] = rng.randrange(1 << dimension)
        else:
            coordinate = rng.randrange(dimension)
            bit = 1 << coordinate
            ones = [index for index, row in enumerate(current) if row & bit]
            zeros = [index for index, row in enumerate(current) if not row & bit]
            if not ones or not zeros:
                continue
            current[rng.choice(ones)] ^= bit
            current[rng.choice(zeros)] |= bit
        if not admissible_rows(current):
            current[:] = previous
            continue

        key = tuple(sorted(current))
        candidate_score, candidate_a, candidate_bad = score_of(key)
        evaluated += 1
        if candidate_a < 0.0:
            recorded = negatives.get(key)
            if recorded is None or candidate_bad < recorded[0]:
                negatives[key] = (candidate_bad, candidate_a)
            if candidate_bad < best[0]:
                best = (candidate_bad, key, candidate_a)
        delta = candidate_score - current_score
        if delta <= 0 or rng.random() < 2.0 ** (-delta / temperature):
            current_score, current_a, current_bad = (
                candidate_score, candidate_a, candidate_bad)
            accepted += 1
        else:
            current[:] = previous

    # Re-evaluate on all n! orders only the most promising witnesses: the walk
    # can visit thousands of negatives, and a full n=7 average costs 0.74 s.
    # Ranking by defect first keeps exactly the ones the frontier needs.
    ranked = sorted(negatives.items(), key=lambda kv: kv[1][0])[:60]
    frontier = []
    for rows, (bad, screened) in ranked:
        full_a, _q, _c = full_objective(rows, dimension)
        if full_a >= 0.0:
            continue
        exact_defect = Fraction(bad, denominator)
        assert admissible(rows, dimension)
        frontier.append({
            "family": list(rows),
            "defect": str(exact_defect),
            "defect_decimal": float(exact_defect),
            "A_plus": full_a,
            "screened_A_plus": screened,
            "ratio": ratio(full_a, exact_defect),
            "normalized": normalized(rows, dimension),
        })
    frontier.sort(key=lambda row: row["defect_decimal"])
    return {
        "seeded": True,
        "size": size,
        "steps": steps,
        "margin": margin,
        "penalty": penalty,
        "objective_orders": "all" if use_full else screen,
        "evaluated": evaluated,
        "accepted": accepted,
        "start_defect": str(Fraction(defect_count(start), denominator)),
        "negatives_visited": len(frontier),
        "lowest_defect_negative": frontier[0] if frontier else None,
        "best_ratio_negative": (
            max(frontier, key=lambda row: row["ratio"]) if frontier else None
        ),
        "frontier": frontier[:40],
    }


def hunt(
    dimension: int,
    size: int,
    lam: float,
    rng: Random,
    steps: int,
    screen: int,
) -> dict[str, object]:
    orders = tuple(
        tuple(rng.sample(range(dimension), dimension)) for _ in range(screen)
    )
    degrees = degree_sequences(size, dimension, rng)
    rows = seed_family(size, dimension, degrees, rng)
    if rows is None:
        return {"seeded": False}

    cap = 2 * size // 5
    need = reimer_threshold(size)
    denominator = size ** 2

    def admissible_rows(candidate: list[int]) -> bool:
        if len(set(candidate)) != len(candidate):
            return False
        counts = coordinate_counts(candidate, dimension)
        return max(counts) <= cap and sum(counts) >= need

    def score_of(candidate: list[int]) -> tuple[float, float, int]:
        key = tuple(sorted(candidate))
        bad = defect_count(key)
        a_plus, _iid, _one = objective(key, dimension, orders)
        return a_plus + lam * bad / denominator, a_plus, bad

    current = list(rows)
    current_score, current_a, current_bad = score_of(current)
    best = (current_score, tuple(sorted(current)), current_a, current_bad)
    evaluated = 1
    accepted = 0
    hot = max(1e-3, abs(lam) / 8.0 + 0.01)

    for step in range(steps):
        temperature = hot * (1.0 - step / steps) + 1e-9
        previous = list(current)
        draw = rng.random()

        if draw < 0.15:
            members = set(current)
            holes = [
                left | right
                for left in current for right in current
                if (left | right) not in members
            ]
            if not holes:
                break
            current[rng.randrange(size)] = rng.choice(holes)
        elif draw < 0.30:
            current[rng.randrange(size)] = rng.randrange(1 << dimension)
        else:
            coordinate = rng.randrange(dimension)
            bit = 1 << coordinate
            ones = [index for index, row in enumerate(current) if row & bit]
            zeros = [index for index, row in enumerate(current) if not row & bit]
            if not ones or not zeros:
                continue
            current[rng.choice(ones)] ^= bit
            current[rng.choice(zeros)] |= bit

        if not admissible_rows(current):
            current[:] = previous
            continue

        candidate_score, candidate_a, candidate_bad = score_of(current)
        evaluated += 1
        delta = candidate_score - current_score
        if delta <= 0 or rng.random() < 2.0 ** (-delta / temperature):
            current_score, current_a, current_bad = (
                candidate_score, candidate_a, candidate_bad)
            accepted += 1
            if candidate_score < best[0]:
                best = (candidate_score, tuple(sorted(current)),
                        candidate_a, candidate_bad)
        else:
            current[:] = previous

    score, family, screen_a, bad = best
    assert admissible(family, dimension), (dimension, family)
    assert defect_count(family) == bad
    exact_defect = Fraction(bad, denominator)
    full_a, full_q, full_c = full_objective(family, dimension)
    return {
        "seeded": True,
        "size": size,
        "degrees": list(degrees),
        "steps": steps,
        "screen_orders": screen,
        "evaluated": evaluated,
        "accepted": accepted,
        "family": list(family),
        "normalized": normalized(family, dimension),
        "defect_count": bad,
        "defect": str(exact_defect),
        "defect_decimal": float(exact_defect),
        "screen_score": score,
        "screen_A_plus": screen_a,
        "A_plus": full_a,
        "Q": full_q,
        "C_plus": full_c,
        "full_score": full_a + lam * float(exact_defect),
        "ratio": ratio(full_a, exact_defect),
        "negative": full_a < 0.0,
    }


def read_checkpoint(path: Path):
    done: dict[tuple[int, str, int, int], dict] = {}
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
            key = (
                payload["dimension"], payload["lam"],
                payload["size"], payload["restart"],
            )
            done[key] = payload
    return done, foreign, truncated


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, default=7)
    parser.add_argument("--sizes", type=str, default="")
    parser.add_argument("--lam", type=str, default="0.04,0.06,0.10")
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument("--screen", type=int, default=SCREEN_ORDERS)
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument(
        "--mode", choices=("hunt", "descend", "sweep"), default="hunt")
    parser.add_argument(
        "--defect-caps", type=str, default="0.70,0.60,0.50,0.45,0.40,0.35")
    parser.add_argument(
        "--defect-checkpoint", type=Path,
        default=HERE / "experiments" / "local_defect_checkpoint.jsonl")
    parser.add_argument("--margin", type=float, default=1e-4)
    parser.add_argument("--penalty", type=float, default=50.0)
    args = parser.parse_args()

    dimension = args.dimension
    if args.sizes:
        sizes = tuple(int(token) for token in args.sizes.split(","))
    else:
        top = max(
            size for size in range(3, (1 << dimension) + 1)
            if reimer_threshold(size) <= dimension * (2 * size // 5)
        )
        sizes = (top,)
    lams = tuple(float(token) for token in args.lam.split(","))

    done, foreign, truncated = read_checkpoint(args.checkpoint)
    started = time.monotonic()
    resumed = 0
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "sweep":
        from search_local_defect import read_witnesses

        witnesses = read_witnesses(args.defect_checkpoint, dimension)
        sources = (args.defect_checkpoint, args.checkpoint)
        rows = []
        for size in sizes:
            witness = witnesses.get(size)
            if witness is None:
                print("no witness for m=%d" % size, file=sys.stderr)
                continue
            pool = seed_pool(sources, dimension, size)
            try:
                base = known_base(dimension)
            except ValueError:
                base = None
            if base is not None and len(base) == size:
                pool.append((base, defect_count(base), None))
            floor = Fraction(witness["defect_count"], size ** 2)
            for token in args.defect_caps.split(","):
                cap_value = Fraction(token)
                if cap_value < floor:
                    continue
                seed_rows = best_seed(pool, cap_value, size)
                if seed_rows is None:
                    seed_rows = tuple(witness["family"])
                for restart in range(args.restarts):
                    key = (dimension, "sweep2:" + token, size, restart)
                    if key in done:
                        rows.append(done[key])
                        resumed += 1
                        continue
                    rng = Random(
                        restart_seed(args.seed, dimension, size, restart)
                        + int(float(token) * 1_000_000) * 15_487
                    )
                    outcome = sweep(
                        dimension, size, cap_value, seed_rows,
                        rng, args.steps, args.screen,
                    )
                    payload = {
                        "schema": RECORD_SCHEMA,
                        "dimension": dimension,
                        "lam": "sweep2:" + token,
                        "size": size,
                        "restart": restart,
                        "seed": args.seed,
                        **outcome,
                    }
                    with args.checkpoint.open("a") as handle:
                        handle.write(json.dumps(payload, sort_keys=True) + "\n")
                        handle.flush()
                    done[key] = payload
                    rows.append(payload)
                    print(
                        "n=%d m=%d cap<=%s restart %d: min_A+=%.6f defect=%s%s (%.0f s)"
                        % (dimension, size, token, restart,
                           payload["min_A_plus"], payload["defect"],
                           "  NEGATIVE" if payload["negative"] else "",
                           time.monotonic() - started),
                        file=sys.stderr, flush=True,
                    )

        curve = {}
        for payload in rows:
            token = payload["lam"].split(":", 1)[1]
            entry = curve.setdefault(token, payload)
            if payload["min_A_plus"] < entry["min_A_plus"]:
                curve[token] = payload
        ordered = sorted(curve.values(), key=lambda row: -row["defect_cap_decimal"])
        negatives = [row for row in ordered if row["negative"]]
        result = {
            "mode": "sweep",
            "dimension": dimension,
            "sizes": list(sizes),
            "restarts": args.restarts,
            "steps": args.steps,
            "resumed_restarts": resumed,
            "foreign_schema_lines": foreign,
            "truncated_checkpoint_lines": truncated,
            "runtime_seconds": round(time.monotonic() - started, 2),
            "curve": [
                {
                    "defect_cap": row["defect_cap"],
                    "defect_cap_decimal": row["defect_cap_decimal"],
                    "min_A_plus": row["min_A_plus"],
                    "defect": row["defect"],
                    "ratio": row["ratio"],
                    "negative": row["negative"],
                    "size": row["size"],
                }
                for row in ordered
            ],
            "lowest_negative_defect_cap": (
                min(negatives, key=lambda row: row["defect_cap_decimal"])["defect_cap"]
                if negatives else None
            ),
            "coverage": (
                "SEEDED HEURISTIC upper bound on Lambda_m(e); exact defect; "
                "A_+ recomputed over all n! orders for every reported point"
            ),
        }
        encoded = json.dumps(result, indent=2, sort_keys=True)
        if args.report is not None:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(encoded + "\n")
        print(encoded)
        return

    if args.mode == "descend":
        records = []
        for restart in range(args.restarts):
            key = (dimension, "descend", 0, restart)
            if key in done:
                records.append(done[key])
                resumed += 1
                continue
            rng = Random(restart_seed(args.seed, dimension, 0, restart) + 104_729)
            outcome = descend(
                dimension, rng, args.steps, args.margin,
                args.penalty, args.screen,
            )
            outcome.pop("seeded", None)
            payload = {
                "schema": RECORD_SCHEMA,
                "dimension": dimension,
                "lam": "descend",
                "size": 0,
                "restart": restart,
                "seed": args.seed,
                **outcome,
            }
            with args.checkpoint.open("a") as handle:
                handle.write(json.dumps(payload, sort_keys=True) + "\n")
                handle.flush()
            done[key] = payload
            records.append(payload)
            low = payload.get("lowest_defect_negative")
            print(
                "n=%d descend restart %d/%d: negatives=%d lowest_defect=%s "
                "ratio=%s (%.0f s)"
                % (dimension, restart + 1, args.restarts,
                   payload.get("negatives_visited", 0),
                   None if low is None else low["defect"],
                   None if low is None else round(low["ratio"], 6),
                   time.monotonic() - started),
                file=sys.stderr, flush=True,
            )

        witnesses = [
            row for payload in records for row in payload.get("frontier", [])
        ]
        witnesses.sort(key=lambda row: row["defect_decimal"])
        best_ratio = max(witnesses, key=lambda row: row["ratio"], default=None)
        result = {
            "mode": "descend",
            "dimension": dimension,
            "restarts": args.restarts,
            "steps": args.steps,
            "margin": args.margin,
            "penalty": args.penalty,
            "resumed_restarts": resumed,
            "foreign_schema_lines": foreign,
            "truncated_checkpoint_lines": truncated,
            "runtime_seconds": round(time.monotonic() - started, 2),
            "negative_witnesses": len(witnesses),
            "lowest_defect_negative": witnesses[0] if witnesses else None,
            "best_ratio_negative": best_ratio,
            "published_base_defect": (
                "444/625" if dimension == 6 else "64/81"
            ),
            "coverage": (
                "SEEDED HEURISTIC; exact defect; A_+ float64 over all n! orders "
                "for every reported witness"
            ),
        }
        encoded = json.dumps(result, indent=2, sort_keys=True)
        if args.report is not None:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(encoded + "\n")
        print(encoded)
        return

    for lam in lams:
        for size in sizes:
            for restart in range(args.restarts):
                key = (dimension, str(lam), size, restart)
                if key in done:
                    resumed += 1
                    continue
                rng = Random(
                    restart_seed(args.seed, dimension, size, restart)
                    + int(lam * 1_000_000) * 7_919
                )
                outcome = hunt(
                    dimension, size, lam, rng, args.steps, args.screen)
                if not outcome.pop("seeded"):
                    outcome = {"size": size, "ratio": None}
                payload = {
                    "schema": RECORD_SCHEMA,
                    "dimension": dimension,
                    "lam": str(lam),
                    "size": size,
                    "restart": restart,
                    "seed": args.seed,
                    **outcome,
                }
                with args.checkpoint.open("a") as handle:
                    handle.write(json.dumps(payload, sort_keys=True) + "\n")
                    handle.flush()
                done[key] = payload
                print(
                    "n=%d m=%d lam=%s restart %d/%d: A+=%s defect=%s ratio=%s (%.0f s)"
                    % (dimension, size, lam, restart + 1, args.restarts,
                       None if outcome.get("A_plus") is None
                       else round(outcome["A_plus"], 6),
                       outcome.get("defect"),
                       None if outcome.get("ratio") is None
                       else round(outcome["ratio"], 6),
                       time.monotonic() - started),
                    file=sys.stderr, flush=True,
                )

    best_ratio = None
    best_negative_defect = None
    rows = []
    for (dim, lam, size, restart), payload in sorted(done.items()):
        if dim != dimension:
            continue
        rows.append(payload)
        value = payload.get("ratio")
        if value is not None and (best_ratio is None or value > best_ratio["ratio"]):
            best_ratio = payload
        if payload.get("negative") and (
            best_negative_defect is None
            or payload["defect_decimal"] < best_negative_defect["defect_decimal"]
        ):
            best_negative_defect = payload

    result = {
        "mode": "hunt",
        "dimension": dimension,
        "lams": list(lams),
        "sizes": list(sizes),
        "restarts": args.restarts,
        "steps": args.steps,
        "screen_orders": args.screen,
        "resumed_restarts": resumed,
        "foreign_schema_lines": foreign,
        "truncated_checkpoint_lines": truncated,
        "runtime_seconds": round(time.monotonic() - started, 2),
        "evaluated_records": len(rows),
        "negatives_found": sum(1 for row in rows if row.get("negative")),
        "best_ratio": best_ratio,
        "lowest_defect_negative": best_negative_defect,
        "calibration": {
            "B_n7_ratio": 0.0362591270,
            "D_n6_ratio": 0.0192456471,
            "note": "published Gate B bases; beat B to strengthen the finite base",
        },
        "coverage": "SEEDED HEURISTIC; float64 A_+; exact defect",
    }
    encoded = json.dumps(result, indent=2, sort_keys=True)
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded + "\n")
    print(encoded)


if __name__ == "__main__":
    main()
