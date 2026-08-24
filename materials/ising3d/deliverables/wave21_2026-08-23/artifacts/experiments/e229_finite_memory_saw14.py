#!/usr/bin/env python3
"""Bounded compact crosswalk for the finite-memory SAW automata, k <= 12.

This producer does not construct the memory-14 automaton and makes no
memory-14 endpoint claim.  Its purpose is narrower: rebuild the e227 controls
with packed signed-axis words, prove the representation crosswalk on every
reachable state and transition, and emit a small durable artifact.

The old draft first enumerated every simple suffix through length nine and
then applied all 48 cubic symmetries to each suffix.  That exhaustive preamble
is unnecessary.  A signed-permutation orbit has a unique *first-use normal
form*: axes are renamed 0, 1, 2 in order of first occurrence and the first
step on each new axis is made positive.  e227's lexicographically least
coordinate representative is exactly the coordinatewise negative of this
normal form.  The producer applies that bijection lazily only to states reached
by the k <= 12 automata, then matches e227's stored state and transition
SHA-256 digests exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import resource
import sys
import time
from collections import defaultdict, deque
from datetime import datetime, timezone
from decimal import Decimal, ROUND_CEILING
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "bounds" / "finite_memory_saw14.json"
CONTROL_ARTIFACT = ROOT / "results" / "bounds" / "finite_memory_saw.json"
MEMORIES = (4, 6, 8, 10, 12)
DIMENSION = 3
DIRECTIONS = tuple(
    tuple(sign if coordinate == axis else 0 for coordinate in range(DIMENSION))
    for axis in range(DIMENSION)
    for sign in (1, -1)
)
# e227 scans +axis before -axis.  The crosswalk negates compact coordinates,
# so this reversed signed-axis scan preserves e227's BFS discovery order.
EXTENSION_STEPS = (1, 0, 3, 2, 5, 4)
EXPECTED_COUNTS = {
    4: (3, 7, 14),
    6: (20, 69, 88),
    8: (205, 805, 871),
    10: (2722, 11074, 11365),
    12: (41424, 169975, 171448),
}
PUBLISHED_RHO_CEILINGS = {
    4: Decimal("4.8646"),
    6: Decimal("4.8075"),
    8: Decimal("4.7780"),
    10: Decimal("4.7599"),
    12: Decimal("4.7476"),
}
DEFAULT_CPU_BUDGET_SECONDS = 300.0
RSS_CAP_BYTES = 2_000_000_000
CLOSABILITY_CACHE_LIMIT = 262_144

StepWord = tuple[int, ...]
Point = tuple[int, int, int]
State = tuple[Point, ...]
Row = tuple[tuple[int, int], ...]


class ResourceWall(RuntimeError):
    """A checked process-CPU or resident-memory budget was reached."""

    def __init__(self, kind: str, detail: dict[str, object]):
        self.kind = kind
        self.detail = detail
        super().__init__(f"{kind}: {detail}")


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def resource_snapshot(started: float) -> dict[str, object]:
    return {
        "process_cpu_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }


def budget_tick(started: float, cpu_budget_seconds: float, stage: str) -> None:
    usage = resource_snapshot(started)
    if float(usage["process_cpu_seconds"]) > cpu_budget_seconds:
        raise ResourceWall("process_cpu_budget", {"stage": stage, **usage})
    if int(usage["peak_rss_bytes"]) >= RSS_CAP_BYTES:
        raise ResourceWall("rss_cap", {"stage": stage, **usage})


# ---------------------------------------------------------------- packed words
def pack(steps: Iterable[int]) -> int:
    """Pack at most fifteen 3-bit signed-axis steps plus a 4-bit length."""

    word = 0
    length = 0
    for step in steps:
        if not 0 <= step < 6:
            raise ValueError(f"invalid signed-axis step {step}")
        if length == 15:
            raise ValueError("compact word is limited to fifteen steps")
        word |= int(step) << (3 * length + 5)
        length += 1
    return word | (length << 1) | 1


def length_of(word: int) -> int:
    return (word >> 1) & 0xF


def step_at(word: int, index: int) -> int:
    return (word >> (3 * index + 5)) & 0x7


def steps_of(word: int) -> StepWord:
    return tuple(step_at(word, index) for index in range(length_of(word)))


def positions_of(word: int) -> State:
    point = [0, 0, 0]
    points: list[Point] = [(0, 0, 0)]
    for step in steps_of(word):
        point[step >> 1] += -1 if step & 1 else 1
        points.append((point[0], point[1], point[2]))
    return tuple(points)


def canonical_word(steps: StepWord) -> int:
    """Return the unique first-use normal form of a signed-axis word.

    The first old axis encountered becomes output axis 0, the next unseen old
    axis becomes output axis 1, and so on.  Each axis sign is chosen so its
    first occurrence is positive.  Axis permutations and independent sign
    flips therefore have one and only one image in this normal form.
    """

    axis_map = [-1, -1, -1]
    first_sign = [0, 0, 0]
    next_axis = 0
    normalized: list[int] = []
    for step in steps:
        old_axis = step >> 1
        negative = step & 1
        new_axis = axis_map[old_axis]
        if new_axis < 0:
            new_axis = next_axis
            next_axis += 1
            axis_map[old_axis] = new_axis
            first_sign[old_axis] = negative
        normalized.append(2 * new_axis + (negative ^ first_sign[old_axis]))
    return pack(normalized)


def is_first_use_normal(word: int) -> bool:
    next_axis = 0
    seen: set[int] = set()
    for step in steps_of(word):
        axis = step >> 1
        if axis not in seen:
            if axis != next_axis or step & 1:
                return False
            seen.add(axis)
            next_axis += 1
    return True


def e227_state_from_word(word: int) -> State:
    """Map a normal word to e227's lexicographically least coordinate state.

    Before a source axis is first used, every image on all still-free output
    axes is zero.  At its first use, lexicographic minimization must put it on
    the earliest free output axis and make that first coordinate negative.
    Thus e227's representative is exactly the global coordinate negation of
    the first-use normal form; no 48-way enumeration is needed.
    """

    if not is_first_use_normal(word):
        raise AssertionError("crosswalk requires a first-use-normal word")
    point = [0, 0, 0]
    state: list[Point] = [(0, 0, 0)]
    for step in steps_of(word):
        point[step >> 1] += 1 if step & 1 else -1
        state.append((point[0], point[1], point[2]))
    return tuple(state)


# --------------------------------------------------------------- closability
def append_step(word: int, tail: tuple[Point, ...], step: int) -> tuple[StepWord, tuple[Point, ...]] | None:
    endpoint = tail[-1]
    axis = step >> 1
    candidate = list(endpoint)
    candidate[axis] += -1 if step & 1 else 1
    point = (candidate[0], candidate[1], candidate[2])
    if point == (0, 0, 0) or point in tail:
        return None
    return steps_of(word) + (step,), tail + (point,)


def path_closable(points: State, memory: int) -> bool:
    """Exact finite BFS for a simple completion of total length <= memory."""

    budget = memory - (len(points) - 1)
    if budget < 1:
        return False
    start = points[0]
    endpoint = points[-1]
    blocked = frozenset(points[1:-1])
    queue: deque[tuple[Point, int]] = deque([(endpoint, 0)])
    seen = {endpoint}
    while queue:
        point, distance = queue.popleft()
        if distance >= budget:
            continue
        for direction in DIRECTIONS:
            candidate = (
                point[0] + direction[0],
                point[1] + direction[1],
                point[2] + direction[2],
            )
            if candidate == start:
                return True
            if candidate in blocked or candidate in seen:
                continue
            next_distance = distance + 1
            remaining = sum(abs(candidate[axis] - start[axis]) for axis in range(DIMENSION))
            if next_distance + remaining > budget:
                continue
            seen.add(candidate)
            queue.append((candidate, next_distance))
    return False


class ClosabilityCache:
    """A deterministic capped cache of exact verdicts for normal suffixes."""

    def __init__(self, limit: int = CLOSABILITY_CACHE_LIMIT) -> None:
        self.limit = limit
        self.values: dict[int, bool] = {}
        self.lookups = 0
        self.hits = 0
        self.inserted = 0
        self.uncached_misses = 0
        self.max_suffix_length = 0

    def closable(self, word: int, memory: int) -> bool:
        self.lookups += 1
        self.max_suffix_length = max(self.max_suffix_length, length_of(word))
        value = self.values.get(word)
        if value is not None:
            self.hits += 1
            return value
        value = path_closable(positions_of(word), memory)
        if len(self.values) < self.limit:
            self.values[word] = value
            self.inserted += 1
        else:
            self.uncached_misses += 1
        return value

    def stats(self) -> dict[str, int]:
        return {
            "limit": self.limit,
            "resident_entries": len(self.values),
            "lookups": self.lookups,
            "hits": self.hits,
            "misses": self.lookups - self.hits,
            "inserted": self.inserted,
            "uncached_misses_after_full": self.uncached_misses,
            "max_suffix_length": self.max_suffix_length,
        }


def reduce_state(extension: StepWord, memory: int, cache: ClosabilityCache) -> int:
    for dropped in range(len(extension) + 1):
        suffix = canonical_word(extension[dropped:])
        if cache.closable(suffix, memory):
            return suffix
    raise AssertionError("the empty suffix is always closable")


# ------------------------------------------------------------------ automaton
def build_automaton(
    memory: int, *, run_started: float, cpu_budget_seconds: float
) -> tuple[list[int], list[Row], dict[str, int]]:
    if memory not in MEMORIES:
        raise ValueError(f"compact control is scoped to {MEMORIES}, got {memory}")
    cache = ClosabilityCache()
    seed = canonical_word((0,))
    states = [seed]
    tails = [positions_of(seed)[1:]]
    index = {seed: 0}
    rows: list[Row] = []
    queue: deque[int] = deque([0])
    while queue:
        source = queue.popleft()
        state = states[source]
        tail = tails[source]
        targets: dict[int, int] = defaultdict(int)
        for step in EXTENSION_STEPS:
            extended = append_step(state, tail, step)
            if extended is None:
                continue
            extension, _ = extended
            target = reduce_state(extension, memory, cache)
            target_index = index.get(target)
            if target_index is None:
                target_index = len(states)
                index[target] = target_index
                states.append(target)
                tails.append(positions_of(target)[1:])
                queue.append(target_index)
            targets[target_index] += 1
        rows.append(tuple(sorted(targets.items())))
        if len(rows) % 1024 == 0:
            budget_tick(
                run_started,
                cpu_budget_seconds,
                f"memory {memory}: {len(rows)}/{len(states)} rows/states",
            )
    if len(rows) != len(states) or len(tails) != len(states):
        raise AssertionError("state, geometry, and transition table lengths disagree")
    if not all(is_first_use_normal(state) for state in states):
        raise AssertionError("builder admitted a non-normal compact state")
    return states, rows, cache.stats()


# --------------------------------------------------------------------- digests
def sha256_ascii(payload: str) -> str:
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def state_tag(state: State) -> str:
    return "/".join(",".join(str(value) for value in point) for point in state)


def compact_state_digest(states: Iterable[int]) -> str:
    return sha256_ascii(";".join(str(state) for state in states))


def e227_state_digest(states: Iterable[State]) -> str:
    return sha256_ascii(";".join(state_tag(state) for state in states))


def e227_transition_digest(rows: Iterable[Row]) -> str:
    return sha256_ascii(
        ";".join(",".join(f"{target}:{multiplicity}" for target, multiplicity in row) for row in rows)
    )


def state_content_digest(states: Iterable[State]) -> str:
    payload = json.dumps(sorted(states), separators=(",", ":"))
    return sha256_ascii(payload)


def transition_content_digest(states: list[State], rows: list[Row]) -> str:
    records: list[tuple[State, tuple[tuple[State, int], ...]]] = []
    for source, row in enumerate(rows):
        targets = tuple(sorted((states[target], multiplicity) for target, multiplicity in row))
        records.append((states[source], targets))
    payload = json.dumps(sorted(records), separators=(",", ":"))
    return sha256_ascii(payload)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def published_ceiling(value: object) -> Decimal:
    decimal_value = Decimal(str(value))
    return (decimal_value * Decimal(10_000)).to_integral_value(rounding=ROUND_CEILING) / Decimal(10_000)


def append_check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    if any(check["name"] == name for check in checks):
        raise AssertionError(f"duplicate check name: {name}")
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


# ---------------------------------------------------------------------- runner
def run(cpu_budget_seconds: float = DEFAULT_CPU_BUDGET_SECONDS) -> dict[str, object]:
    if not math.isfinite(cpu_budget_seconds) or cpu_budget_seconds <= 0:
        raise ValueError("process CPU budget must be a positive finite number")
    started = time.process_time()
    legacy_payload = json.loads(CONTROL_ARTIFACT.read_text())
    legacy_rows = {
        int(row["memory"]): row for row in legacy_payload["data"]["memories"]
    }
    if tuple(sorted(legacy_rows)) != MEMORIES:
        raise AssertionError(f"e227 control memories changed: {tuple(sorted(legacy_rows))}")

    checks: list[dict[str, object]] = []
    memory_rows: list[dict[str, object]] = []
    for memory in MEMORIES:
        budget_tick(started, cpu_budget_seconds, f"before memory {memory}")
        build_started = time.process_time()
        rss_before = max_rss_bytes()
        states, rows, cache_stats = build_automaton(
            memory,
            run_started=started,
            cpu_budget_seconds=cpu_budget_seconds,
        )
        build_cpu = time.process_time() - build_started
        rss_after = max_rss_bytes()
        tuple_states = [e227_state_from_word(state) for state in states]
        counts = (
            len(states),
            sum(len(row) for row in rows),
            sum(multiplicity for row in rows for _, multiplicity in row),
        )
        expected_counts = EXPECTED_COUNTS[memory]
        legacy = legacy_rows[memory]
        state_sha = e227_state_digest(tuple_states)
        transition_sha = e227_transition_digest(rows)
        compact_sha = compact_state_digest(states)
        state_content_sha = state_content_digest(tuple_states)
        transition_content_sha = transition_content_digest(tuple_states, rows)
        stored_counts = (
            int(legacy["state_count"]),
            int(legacy["transition_count"]),
            int(legacy["weighted_transition_count"]),
        )
        observed_published_ceiling = published_ceiling(legacy["perron_float_approximation"])

        append_check(
            checks,
            f"memory_{memory}_exact_counts_match_control",
            counts == expected_counts == stored_counts,
            f"states/transitions/weighted={counts}",
        )
        append_check(
            checks,
            f"memory_{memory}_e227_state_digest_exact",
            state_sha == str(legacy["state_sha256"]),
            f"compact-crosswalk={state_sha}; e227={legacy['state_sha256']}",
        )
        append_check(
            checks,
            f"memory_{memory}_e227_transition_digest_exact",
            transition_sha == str(legacy["transition_sha256"]),
            f"compact-crosswalk={transition_sha}; e227={legacy['transition_sha256']}",
        )
        append_check(
            checks,
            f"memory_{memory}_published_rho_rounding_control",
            observed_published_ceiling == PUBLISHED_RHO_CEILINGS[memory],
            (
                f"exact matrix matches e227; e227 rho={legacy['perron_float_approximation']} "
                f"ceil4={observed_published_ceiling} published={PUBLISHED_RHO_CEILINGS[memory]}"
            ),
        )
        append_check(
            checks,
            f"memory_{memory}_closability_cache_bounded",
            cache_stats["resident_entries"] <= cache_stats["limit"],
            json.dumps(cache_stats, sort_keys=True),
        )
        append_check(
            checks,
            f"memory_{memory}_rss_below_two_gib",
            rss_after < RSS_CAP_BYTES,
            f"peak_rss_bytes={rss_after}; cap={RSS_CAP_BYTES}",
        )
        budget_tick(started, cpu_budget_seconds, f"after memory {memory}")

        memory_rows.append(
            {
                "memory": memory,
                "frontier_status": "compact_crosswalk_complete",
                "state_count": counts[0],
                "transition_count": counts[1],
                "weighted_transition_count": counts[2],
                "compact_state_sha256": compact_sha,
                "e227_format_state_sha256": state_sha,
                "e227_format_transition_sha256": transition_sha,
                "state_content_sha256": state_content_sha,
                "transition_content_sha256": transition_content_sha,
                "e227_control": {
                    "state_sha256": legacy["state_sha256"],
                    "transition_sha256": legacy["transition_sha256"],
                    "perron_float_approximation": legacy["perron_float_approximation"],
                    "published_rho_ceiling_4dp": format(PUBLISHED_RHO_CEILINGS[memory], ".4f"),
                    "spectral_value_recomputed": False,
                    "spectral_transfer_basis": "exact ordered state and transition table equality",
                },
                "resources": {
                    "build_process_cpu_seconds": round(build_cpu, 6),
                    "peak_rss_before_bytes": rss_before,
                    "peak_rss_after_bytes": rss_after,
                },
                "closability_cache": cache_stats,
                "observed_wall": None,
            }
        )

    memory_14 = {
        "memory": 14,
        "launch_status": "not_launched",
        "reason": (
            "This bounded crosswalk is structurally restricted to k<=12; "
            "the explicit task contract forbids a k=14 launch."
        ),
        "observed_wall": None,
        "endpoint_claim": None,
    }
    append_check(
        checks,
        "memory_14_not_launched_by_bounded_producer",
        memory_14["launch_status"] == "not_launched",
        str(memory_14["reason"]),
    )
    append_check(
        checks,
        "memory_14_endpoint_claim_absent",
        memory_14["endpoint_claim"] is None,
        "no memory-14 spectral or Ising endpoint statement is emitted",
    )
    names_before_uniqueness_check = [str(check["name"]) for check in checks]
    append_check(
        checks,
        "artifact_check_names_unique",
        len(names_before_uniqueness_check) == len(set(names_before_uniqueness_check)),
        f"{len(names_before_uniqueness_check)} prior names are pairwise distinct",
    )
    failures = [str(check["name"]) for check in checks if not bool(check["passed"])]
    if failures:
        raise AssertionError(f"compact crosswalk checks failed: {failures}")
    budget_tick(started, cpu_budget_seconds, "artifact assembly")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e229_finite_memory_saw14.py",
            "interpreter": sys.executable,
            "process_cpu_budget_seconds": cpu_budget_seconds,
            "rss_cap_bytes": RSS_CAP_BYTES,
            "process_cpu_seconds": resource_snapshot(started)["process_cpu_seconds"],
            "peak_rss_bytes": max_rss_bytes(),
            "observed_wall": None,
            "e227_control_artifact": "results/bounds/finite_memory_saw.json",
            "e227_control_artifact_sha256": file_sha256(CONTROL_ARTIFACT),
            "arithmetic": "exact integer states, transitions, counts, and SHA-256 crosswalks",
        },
        "data": {
            "dimension": DIMENSION,
            "memories": memory_rows,
            "memory_14": memory_14,
            "scope": {
                "proved": "compact first-use automata equal e227 automata content-wise for k=4,6,8,10,12",
                "not_claimed": "no memory-14 automaton, spectral radius, Collatz certificate, or K_c endpoint",
            },
            "canonical_crosswalk": {
                "compact_normal_form": (
                    "rename axes by first occurrence and flip each source-axis sign so its first step is positive"
                ),
                "e227_coordinate_map": "negate every compact-normal coordinate",
                "control_strategy": (
                    "apply the algebraic bijection lazily to reachable compact states only; "
                    "match e227's ordered state and transition digests"
                ),
                "exhaustive_length_9_preamble_used": False,
            },
        },
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cpu-budget-seconds",
        type=float,
        default=DEFAULT_CPU_BUDGET_SECONDS,
        help="hard checked process-CPU budget for the complete k<=12 crosswalk",
    )
    args = parser.parse_args()
    payload = run(args.cpu_budget_seconds)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_name(f".{OUTPUT.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    temporary.replace(OUTPUT)
    for row in payload["data"]["memories"]:
        print(
            f"[PASS] k={row['memory']:2d}: states={row['state_count']:5d}, "
            f"transitions={row['transition_count']:6d}, "
            f"weighted={row['weighted_transition_count']:6d}, "
            f"cpu={row['resources']['build_process_cpu_seconds']:.3f}s, "
            f"peak_rss={row['resources']['peak_rss_after_bytes']}"
        )
    print("[NOT LAUNCHED] k=14; no endpoint claim")
    print(f"[PASS] {len(payload['checks'])} uniquely named checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
