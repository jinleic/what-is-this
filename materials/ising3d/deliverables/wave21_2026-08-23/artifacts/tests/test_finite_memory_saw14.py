#!/usr/bin/env python3
"""Independent compact finite-memory SAW crosswalk verifier.

No producer imports.  Rebuilds memory 4..12 from suffix geometry using an
independent first-use canonicalization and a separate DFS closability test,
then compares the resulting ordered state/transition digests with e227.  It
also requires that the durable artifact record memory-14 as neither launched
nor endpoint-certified.
"""

from __future__ import annotations

import itertools
import json
import time
from collections import defaultdict, deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "bounds" / "finite_memory_saw14.json"
CONTROL = ROOT / "results" / "bounds" / "finite_memory_saw.json"
DIRECTIONS = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
SYMMETRIES = tuple(
    (perm, signs)
    for perm in itertools.permutations(range(3))
    for signs in itertools.product((-1, 1), repeat=3)
)
FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def add(point, step):
    return (point[0] + step[0], point[1] + step[1], point[2] + step[2])


def canonical_tuple(path):
    """A 48-way raw canonicalization used as a code-path control."""

    best = None
    for perm, signs in SYMMETRIES:
        image = tuple(tuple(signs[i] * point[perm[i]] for i in range(3)) for point in path)
        if best is None or image < best:
            best = image
    assert best is not None
    return best


def first_use_canonical(path):
    """Independent first-use normalization, distinct from e229."""

    def step_axis(step):
        return next(index for index, value in enumerate(step) if value)
    axis_map: dict[int, int] = {}
    axis_sign: dict[int, int] = {}
    out = [(0, 0, 0)]
    for source, target in zip(path, path[1:]):
        delta = tuple(target[i] - source[i] for i in range(3))
        step = next(step for step in DIRECTIONS if step == delta)
        axis = step_axis(step)
        sign = step[axis]
        if axis not in axis_map:
            new_axis = len(axis_map)
            axis_map[axis] = new_axis
            axis_sign[axis] = 1 if sign > 0 else -1
        mapped_axis = axis_map[axis]
        mapped_sign = sign * axis_sign[axis]
        mapped = tuple(mapped_sign if i == mapped_axis else 0 for i in range(3))
        out.append(add(out[-1], mapped))
    return out


def canonical_path(path):
    """Map the independent first-use form to e227's sign convention."""

    normal = first_use_canonical(path)
    return tuple(tuple(-value for value in point) for point in normal)


def can_close(path, memory):
    """Independent DFS with exact budget pruning."""

    budget = memory - (len(path) - 1)
    if budget < 1:
        return False
    start, endpoint = path[0], path[-1]
    interior = frozenset(path[1:-1])
    stack = [(endpoint, 0)]
    seen = {endpoint}
    while stack:
        point, used = stack.pop()
        if used >= budget:
            continue
        for step in DIRECTIONS:
            candidate = add(point, step)
            if candidate == start:
                return True
            if candidate in interior or candidate in seen:
                continue
            remaining = sum(abs(candidate[i] - start[i]) for i in range(3))
            if used + 1 + remaining > budget:
                continue
            seen.add(candidate)
            stack.append((candidate, used + 1))
    return False


def shrink(path, memory):
    for drop in range(len(path)):
        suffix = path[drop:]
        if can_close(suffix, memory):
            return canonical_path(suffix)
    raise AssertionError("the empty suffix is closable")


def build(memory):
    seed = canonical_path(((0, 0, 0), (1, 0, 0)))
    states = [seed]
    index = {seed: 0}
    rows: list[dict[int, int]] = []
    queue = deque([seed])
    while queue:
        state = queue.popleft()
        row: dict[int, int] = defaultdict(int)
        for step in DIRECTIONS:
            endpoint = add(state[-1], step)
            if endpoint in state:
                continue
            target = shrink(state + (endpoint,), memory)
            slot = index.get(target)
            if slot is None:
                slot = len(states)
                index[target] = slot
                states.append(target)
                queue.append(target)
            row[slot] += 1
        rows.append(dict(row))
    return states, rows


def tuple_times(states, i, j):
    return tuple(states[j][k] - states[i][k] for k in range(3))


def digest_states(states):
    import hashlib

    return hashlib.sha256(
        ";".join("/".join(",".join(str(c) for c in point) for point in state) for state in states).encode()
    ).hexdigest()


def digest_rows(states, rows):
    import hashlib

    return hashlib.sha256(
        ";".join(",".join(f"{target}:{value}" for target, value in sorted(row.items())) for row in rows).encode()
    ).hexdigest()


def main() -> int:
    began = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    legacy = {int(row["memory"]): row for row in json.loads(CONTROL.read_text())["data"]["memories"]}
    check(
        "artifact recorded memory-14 as not launched with no endpoint",
        artifact["data"]["memory_14"]["launch_status"] == "not_launched"
        and artifact["data"]["memory_14"]["endpoint_claim"] is None,
    )
    check(
        "artifact checks uniquely named and all true",
        len({record["name"] for record in artifact["checks"]}) == len(artifact["checks"])
        and all(record["passed"] for record in artifact["checks"]),
        f"{len(artifact['checks'])} checks",
    )
    for memory in (4, 6, 8, 10, 12):
        states, rows = build(memory)
        stored = legacy[memory]
        check(
            f"memory_{memory}_independent_compact_rebuild_matches_content",
            len(states) == int(stored["state_count"])
            and sum(len(row) for row in rows) == int(stored["transition_count"])
            and sum(sum(row.values()) for row in rows) == int(stored["weighted_transition_count"]),
            (
                f"states={len(states)}/transitions={sum(len(row) for row in rows)}/"
                f"weighted={sum(sum(row.values()) for row in rows)}; "
                f"state={'yes' if digest_states(states)==stored['state_sha256'] else 'no'}; "
                f"rows={'yes' if digest_rows(states, rows)==stored['transition_sha256'] else 'no'}"
            ),
        )
    # Raw 48-way control, intentionally different from first-use normalization.
    for memory in (4, 6):
        states, _ = build(memory)
        check(
            f"memory_{memory}_raw_48_way_canonicalization_agrees",
            all(canonical_tuple(state) == state for state in states),
            f"{len(states)} states",
        )
    check("process_budget", time.process_time() - began < 900, f"CPU={time.process_time()-began:.3f}")
    if FAILS:
        print(f"FAIL: {len(FAILS)} checks")
        return 1
    print("PASS: independent compact SAW crosswalk verifier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
