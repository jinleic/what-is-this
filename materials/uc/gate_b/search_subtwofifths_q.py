#!/usr/bin/env python3
"""Search the sub-2/5 defect region for Gate B and clone-reachable witnesses.

The walk always keeps

    eps_vee < 2/5,  [n] in F,  cap, Reimer, distinct rows

as hard constraints.  Its three discovery scores serve different purposes:

* ``--score q`` minimizes sampled iid ``Q``; full-order ``Q`` is then enclosed
  exactly as a necessary rejection gate because ``C_+ >= 0``.
* ``--score a`` minimizes the sampled full order-average objective.
* ``--score min-a`` minimizes the least sampled fixed-order objective.  A
  selected order is then enclosed exactly.  By the sharp clone-order law in
  ``audit_clone_limit.py``, one certified negative fixed order is sufficient
  for an explicit finite clone multiset with the same rows and defect.

The observed sub-2/5 Q survivor shows why ``q`` and either Bellman-aware score
must not be confused: its Q-only room is ``-0.00307`` but sampled
``C_+ ~= 4.82`` contributes ``+0.17168``.  Full-order Q remains a cheap
necessary gate, not the main ranking objective.

Float64 and sampled-order comparisons are discovery only.  The full-order Q
enclosure and the selected fixed-order Bellman enclosure use exact rational
two-sided arithmetic.  Checkpoints are append-only, restart-indexed,
deterministic, and contain no timing or resume-dependent fields.

Run the Bellman-aware modes from an existing defect or Q checkpoint:

    OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
      uc/gate_b/search_subtwofifths_q.py --dimension 9 --sizes 15,20,30 \
      --score min-a --defect-checkpoint ... --extra-seed-checkpoint ...
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb, log2
from pathlib import Path
from random import Random
import sys
from typing import Callable

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from search_local_defect import (  # noqa: E402
    admissible,
    coordinate_counts,
    defect_count,
    pool_families,
    normalized,
    read_witnesses,
    reimer_threshold,
    restart_seed,
)
from verify_gate_b_rational import (  # noqa: E402
    ACCUMULATOR_BITS,
    ALPHA,
    _predecessor_law,
    bellman_bounds,
    decimal_lower,
    decimal_upper,
    log2_lower,
    log2_upper,
    round_down,
    round_up,
    shapley_iid_bounds,
)
from audit_clone_limit import (  # noqa: E402
    finite_negative_clone_witness,
    fixed_order_a_bounds,
)

SCHEMA = 1
DEFAULT_CHECKPOINT = HERE / "experiments" / "subtwofifths_q_checkpoint.jsonl"
DEFAULT_REPORT = HERE / "candidates" / "subtwofifths_q_frontier.json"


def entropy(value: Fraction) -> float:
    point = float(value)
    if point <= 0.0 or point >= 1.0:
        return 0.0
    return -point * log2(point) - (1.0 - point) * log2(1.0 - point)


def fixed_order_q_float(
    family: tuple[int, ...], order: tuple[int, ...]
) -> float:
    total = 0.0
    predecessors: list[int] = []
    for coordinate in order:
        law = _predecessor_law(family, coordinate, tuple(predecessors))
        total += sum(
            float(left_mass * right_mass)
            * entropy(left + right - left * right)
            for left, left_mass in law
            for right, right_mass in law
        )
        predecessors.append(coordinate)
    return total


def screen_q(
    family: tuple[int, ...], orders: tuple[tuple[int, ...], ...]
) -> float:
    return sum(fixed_order_q_float(family, order) for order in orders) / len(orders)

def screen_min_a(
    family: tuple[int, ...],
    orders: tuple[tuple[int, ...], ...],
    one_sided_costs: Callable[..., tuple[list[float], object]],
) -> tuple[float, tuple[int, ...]]:
    """Least sampled fixed-order objective and its order; discovery only."""
    q_values = tuple(fixed_order_q_float(family, order) for order in orders)
    c_values, _states = one_sided_costs(
        family, orders, dimension=len(orders[0])
    )
    values = tuple(
        (1.0 - float(ALPHA)) * q_value
        + float(ALPHA) * c_value
        - log2(len(family))
        for q_value, c_value in zip(q_values, c_values)
    )
    index = min(range(len(orders)), key=lambda item: (values[item], orders[item]))
    return values[index], orders[index]


def exact_finalist(
    family: tuple[int, ...],
    dimension: int,
    orders: tuple[tuple[int, ...], ...],
    probe_order: tuple[int, ...] | None = None,
) -> dict[str, object]:
    q_low, q_high = shapley_iid_bounds(family, dimension)
    values = [bellman_bounds(family, dimension, order) for order in orders]
    c_low = round_down(
        sum(value[0] for value in values) / len(values), ACCUMULATOR_BITS
    )
    c_high = round_up(
        sum(value[1] for value in values) / len(values), ACCUMULATOR_BITS
    )
    size = Fraction(len(family))
    q_only_low = round_down(
        (1 - ALPHA) * q_low - log2_upper(size), ACCUMULATOR_BITS
    )
    q_only_high = round_up(
        (1 - ALPHA) * q_high - log2_lower(size), ACCUMULATOR_BITS
    )
    sampled_a_low = round_down(
        (1 - ALPHA) * q_low + ALPHA * c_low - log2_upper(size),
        ACCUMULATOR_BITS,
    )
    sampled_a_high = round_up(
        (1 - ALPHA) * q_high + ALPHA * c_high - log2_lower(size),
        ACCUMULATOR_BITS,
    )
    result = {
        "q_lower": str(q_low),
        "q_upper": str(q_high),
        "q_lower_decimal": decimal_lower(q_low, 12),
        "q_upper_decimal": decimal_upper(q_high, 12),
        "q_only_a_lower": str(q_only_low),
        "q_only_a_upper": str(q_only_high),
        "q_only_a_lower_decimal": decimal_lower(q_only_low, 12),
        "q_only_a_upper_decimal": decimal_upper(q_only_high, 12),
        "q_alone_certifies_positive": q_only_low > 0,
        "sampled_bellman_lower": str(c_low),
        "sampled_bellman_upper": str(c_high),
        "sampled_a_lower": str(sampled_a_low),
        "sampled_a_upper": str(sampled_a_high),
        "sampled_a_lower_decimal": decimal_lower(sampled_a_low, 12),
        "sampled_a_upper_decimal": decimal_upper(sampled_a_high, 12),
        "sampled_a_status": "DISCOVERY ONLY: exact arithmetic over sampled orders",
    }
    if probe_order is not None:
        fixed_low, fixed_high = fixed_order_a_bounds(
            family, dimension, probe_order
        )
        global_lower = round_down(-log2_upper(size), ACCUMULATOR_BITS)
        global_upper = round_up(
            Fraction(dimension) - log2_lower(size), ACCUMULATOR_BITS
        )
        result.update({
            "selected_fixed_order": list(probe_order),
            "selected_fixed_order_a_lower": str(fixed_low),
            "selected_fixed_order_a_upper": str(fixed_high),
            "selected_fixed_order_a_lower_decimal": decimal_lower(fixed_low, 12),
            "selected_fixed_order_a_upper_decimal": decimal_upper(fixed_high, 12),
            "selected_fixed_order_certified_negative": fixed_high < 0,
            "selected_fixed_order_status": (
                "EXACT RATIONAL TWO-SIDED ENCLOSURE; order selected by discovery score"
            ),
            "finite_negative_clone_witness": finite_negative_clone_witness(
                probe_order,
                fixed_low,
                fixed_high,
                global_lower,
                global_upper,
            ),
            "clone_witness_global_order_bounds": {
                "lower": str(global_lower),
                "upper": str(global_upper),
                "derivation": (
                    "-log2(m) <= A_{+,pi} <= n-log2(m) since every entropy "
                    "sum lies in [0,n]"
                ),
            },
        })
    return result

def refresh_clone_witness(payload: dict[str, object]) -> dict[str, object]:
    """Upgrade an append-only finalist record to the current exact clone certificate."""
    order = payload.get("selected_fixed_order")
    if order is None:
        return payload
    refreshed = dict(payload)
    size = Fraction(int(payload["size"]))
    dimension = int(payload["dimension"])
    global_lower = round_down(-log2_upper(size), ACCUMULATOR_BITS)
    global_upper = round_up(
        Fraction(dimension) - log2_lower(size), ACCUMULATOR_BITS
    )
    refreshed["finite_negative_clone_witness"] = finite_negative_clone_witness(
        tuple(order),
        Fraction(payload["selected_fixed_order_a_lower"]),
        Fraction(payload["selected_fixed_order_a_upper"]),
        global_lower,
        global_upper,
    )
    refreshed["clone_witness_global_order_bounds"] = {
        "lower": str(global_lower),
        "upper": str(global_upper),
        "derivation": (
            "-log2(m) <= A_{+,pi} <= n-log2(m) since every entropy "
            "sum lies in [0,n]"
        ),
    }
    return refreshed


def valid(
    rows: list[int], dimension: int, size: int, bad_limit: int
) -> bool:
    if len(set(rows)) != size:
        return False
    if (1 << dimension) - 1 not in rows:
        return False
    counts = coordinate_counts(rows, dimension)
    if max(counts) > 2 * size // 5:
        return False
    if sum(counts) < reimer_threshold(size):
        return False
    return defect_count(tuple(rows)) <= bad_limit


def seed_pool(
    checkpoint: Path,
    dimension: int,
    size: int,
    bad_limit: int,
    extra_checkpoint: Path | None = None,
) -> list[tuple[int, ...]]:
    seen: set[tuple[int, ...]] = set()
    best = read_witnesses(checkpoint, dimension).get(size)
    if best is not None:
        rows = tuple(best["family"])
        if defect_count(rows) <= bad_limit and (1 << dimension) - 1 in rows:
            seen.add(rows)
    for entry in pool_families(checkpoint, dimension, 1000):
        rows = tuple(entry["family"])
        if len(rows) != size:
            continue
        if entry["defect_count"] <= bad_limit and (1 << dimension) - 1 in rows:
            seen.add(rows)
    if extra_checkpoint is not None and extra_checkpoint.exists():
        with extra_checkpoint.open() as handle:
            for line in handle:
                try:
                    payload = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if payload.get("dimension") != dimension or payload.get("size") != size:
                    continue
                rows = tuple(payload.get("family", ()))
                if (
                    len(rows) == size
                    and defect_count(rows) <= bad_limit
                    and (1 << dimension) - 1 in rows
                    and admissible(rows, dimension)
                ):
                    seen.add(rows)
    return sorted(seen)

def strict_bad_limit(size: int, defect_cap: Fraction) -> int:
    """Largest integer failure count strictly below ``defect_cap * size^2``."""
    assert 0 < defect_cap <= 1
    return (
        defect_cap.numerator * size * size - 1
    ) // defect_cap.denominator


def walk(
    start: tuple[int, ...],
    dimension: int,
    scorer: Callable[[tuple[int, ...]], float],
    rng: Random,
    steps: int,
    bad_limit: int,
) -> tuple[tuple[int, ...], float, int, int]:
    size = len(start)
    assert defect_count(start) <= bad_limit
    full = (1 << dimension) - 1
    current = list(start)
    current_score = scorer(tuple(sorted(current)))
    best_rows = tuple(sorted(current))
    best_score = current_score
    accepted = 0
    hot = 0.04
    for step in range(steps):
        previous = list(current)
        draw = rng.random()
        mutable = [index for index, row in enumerate(current) if row != full]
        if draw < 0.20:
            victim = rng.choice(mutable)
            current[victim] = rng.randrange(1 << dimension)
        else:
            coordinate = rng.randrange(dimension)
            bit = 1 << coordinate
            ones = [index for index in mutable if current[index] & bit]
            zeros = [index for index in mutable if not current[index] & bit]
            if not ones or not zeros:
                continue
            source = rng.choice(ones)
            target = rng.choice(zeros)
            current[source] ^= bit
            current[target] |= bit
        if not valid(current, dimension, size, bad_limit):
            current[:] = previous
            continue
        candidate = tuple(sorted(current))
        candidate_score = scorer(candidate)
        temperature = hot * (1.0 - step / steps) + 1e-9
        delta = candidate_score - current_score
        if delta <= 0 or rng.random() < 2.0 ** (-delta / temperature):
            current_score = candidate_score
            accepted += 1
            if candidate_score < best_score:
                best_score = candidate_score
                best_rows = candidate
        else:
            current[:] = previous
    return best_rows, best_score, defect_count(best_rows), accepted


def read_checkpoint(
    path: Path,
) -> dict[tuple[int, int, int, str, int, int, int, str], dict[str, object]]:
    done = {}
    if not path.exists():
        return done
    with path.open() as handle:
        for line in handle:
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            if payload.get("schema") != SCHEMA:
                continue
            required = (
                "dimension", "size", "restart", "score_mode",
                "steps", "screen_order_count", "seed",
            )
            if any(field not in payload for field in required):
                continue
            key = (
                int(payload["dimension"]),
                int(payload["size"]),
                int(payload["restart"]),
                str(payload["score_mode"]),
                int(payload["steps"]),
                int(payload["screen_order_count"]),
                int(payload["seed"]),
                str(payload.get("strict_defect_cap", "2/5")),
            )
            done[key] = payload
    return done

def matching_checkpoint_records(
    done: dict[
        tuple[int, int, int, str, int, int, int, str],
        dict[str, object],
    ],
    dimension: int,
    sizes: tuple[int, ...],
    restarts: int,
    score: str,
    steps: int,
    screen_orders: int,
    seed: int,
    defect_cap: Fraction,
) -> list[dict[str, object]]:
    """Records reproducible by exactly the configuration declared in a report."""
    return [
        refresh_clone_witness(payload)
        for key, payload in sorted(done.items())
        if (
            key[0] == dimension
            and key[1] in sizes
            and key[2] < restarts
            and key[3:] == (
                score,
                steps,
                screen_orders,
                seed,
                str(defect_cap),
            )
        )
    ]


def append(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--sizes", type=str, required=True)
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--steps", type=int, default=1500)
    parser.add_argument("--screen-orders", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--defect-checkpoint", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--score", choices=("q", "a", "min-a"), default="q")
    parser.add_argument(
        "--defect-cap", type=Fraction, default=Fraction(2, 5)
    )
    parser.add_argument("--extra-seed-checkpoint", type=Path, default=None)
    arguments = parser.parse_args()

    dimension = arguments.dimension
    sizes = tuple(int(token) for token in arguments.sizes.split(","))
    done = read_checkpoint(arguments.checkpoint)
    one_sided_costs = None
    if arguments.score in ("a", "min-a"):
        # Delayed: this fast discovery evaluator belongs to the project venv.
        # The exact finalist route above has no scipy/python-flint dependency.
        from shapley_n6_shared_bellman import one_sided_costs as shared_costs
        one_sided_costs = shared_costs
    for size in sizes:
        bad_limit = strict_bad_limit(size, arguments.defect_cap)
        seeds = seed_pool(
            arguments.defect_checkpoint,
            dimension,
            size,
            bad_limit,
            arguments.extra_seed_checkpoint,
        )
        if not seeds:
            print(
                f"n={dimension} m={size}: no dominant seed below "
                f"{arguments.defect_cap}",
                file=sys.stderr,
            )
            continue
        for restart in range(arguments.restarts):
            key = (
                dimension,
                size,
                restart,
                arguments.score,
                arguments.steps,
                arguments.screen_orders,
                arguments.seed,
                str(arguments.defect_cap),
            )
            if key in done:
                continue
            rng = Random(restart_seed(arguments.seed, dimension, size, restart))
            orders = tuple(
                tuple(rng.sample(range(dimension), dimension))
                for _ in range(arguments.screen_orders)
            )
            start = seeds[restart % len(seeds)]

            def scorer(candidate: tuple[int, ...]) -> float:
                q_value = screen_q(candidate, orders)
                if arguments.score == "q":
                    return q_value
                assert one_sided_costs is not None
                if arguments.score == "min-a":
                    return screen_min_a(candidate, orders, one_sided_costs)[0]
                costs, _states = one_sided_costs(
                    candidate, orders, dimension=dimension
                )
                c_value = sum(costs) / len(costs)
                return (
                    (1.0 - float(ALPHA)) * q_value
                    + float(ALPHA) * c_value
                    - log2(len(candidate))
                )

            best, screened_score, bad, accepted = walk(
                start, dimension, scorer, rng, arguments.steps, bad_limit
            )
            assert admissible(best, dimension)
            assert bad <= bad_limit
            probe_order = None
            if arguments.score == "min-a":
                assert one_sided_costs is not None
                _minimum, probe_order = screen_min_a(
                    best, orders, one_sided_costs
                )
            exact = exact_finalist(best, dimension, orders, probe_order)
            payload = {
                "schema": SCHEMA,
                "dimension": dimension,
                "size": size,
                "restart": restart,
                "seed": arguments.seed,
                "steps": arguments.steps,
                "screen_order_count": arguments.screen_orders,
                "accepted": accepted,
                "family": list(best),
                "normalized": normalized(best, dimension),
                "contains_full_set": (1 << dimension) - 1 in best,
                "closure_defect": str(Fraction(bad, size * size)),
                "closure_defect_count": bad,
                "score_mode": arguments.score,
                "strict_defect_cap": str(arguments.defect_cap),
                "screen_score": screened_score,
                **exact,
            }
            append(arguments.checkpoint, payload)
            done[key] = payload
            message = (
                f"n={dimension} m={size} restart={restart} "
                f"defect={payload['closure_defect']} "
                f"Q={payload['q_upper_decimal']} "
                f"Q-only={payload['q_only_a_lower_decimal']} "
                f"sample-A={payload['sampled_a_upper_decimal']}"
            )
            if probe_order is not None:
                message += (
                    " selected-order-A="
                    f"{payload['selected_fixed_order_a_upper_decimal']}"
                )
            print(message, file=sys.stderr, flush=True)

    relevant = matching_checkpoint_records(
        done,
        dimension,
        sizes,
        arguments.restarts,
        arguments.score,
        arguments.steps,
        arguments.screen_orders,
        arguments.seed,
        arguments.defect_cap,
    )

    def ranking_key(payload: dict[str, object]) -> Fraction:
        if arguments.score == "q":
            return Fraction(payload["q_upper"])
        if arguments.score == "a":
            return Fraction(payload["sampled_a_upper"])
        return Fraction(payload["selected_fixed_order_a_upper"])

    by_size = {}
    for size in sizes:
        rows = [payload for payload in relevant if payload["size"] == size]
        if not rows:
            continue
        by_size[str(size)] = min(rows, key=ranking_key)
    frontier = min(by_size.values(), key=ranking_key, default=None)
    negative_fixed_order_finalists = [
        payload
        for payload in relevant
        if payload.get("selected_fixed_order_certified_negative", False)
    ]
    lowest_defect_negative = min(
        negative_fixed_order_finalists,
        key=lambda payload: (
            Fraction(payload["closure_defect"]),
            Fraction(payload["selected_fixed_order_a_upper"]),
        ),
        default=None,
    )
    report = {
        "mode": f"subtwofifths_{arguments.score}_search",
        "coverage": (
            "SEEDED HEURISTIC; exact defect and exact full-order Q enclosure; "
            + (
                "selected fixed-order A_+ has an exact rational two-sided enclosure"
                if arguments.score == "min-a"
                else "Bellman/A_+ sampled-order discovery only"
            )
        ),
        "score_mode": arguments.score,
        "dimension": dimension,
        "sizes": list(sizes),
        "restarts": arguments.restarts,
        "steps": arguments.steps,
        "screen_order_count": arguments.screen_orders,
        "dominant_full_set_required": True,
        "strict_defect_cap": str(arguments.defect_cap),
        "normalization_convention": (
            "active and separating is reported, not required by the displayed supremum"
        ),
        "displayed_supremum_convention": "cap and Reimer only",
        "best_score_finalist": frontier,
        "best_by_size": by_size,
        "lowest_defect_certified_negative_fixed_order": lowest_defect_negative,
        "any_finalist_passes_q_necessary_condition": any(
            not payload["q_alone_certifies_positive"] for payload in relevant
        ),
        "any_finalist_has_certified_negative_fixed_order": any(
            payload.get("selected_fixed_order_certified_negative", False)
            for payload in relevant
        ),
        "clone_reachability_convention": (
            "A certified negative fixed order gives a finite negative clone "
            "multiset in the displayed cap/Reimer class; nontrivial cloning "
            "does not preserve separation"
        ),
    }
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
