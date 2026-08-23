#!/usr/bin/env python3
"""Clean-room verifier for the finite-memory SAW upper bound and K_c corollary.

Independent of experiments/e227_finite_memory_saw.py (imports nothing from it):
  * rebuilds the memory-k automata with a fresh DFS closability test and a
    fresh canonicalization, k = 4..12;
  * checks the exact k=4 characteristic polynomial against the closed-form
    denominator 1 - 4z - 4z^2 - z^3 of Poenitz--Tittmann (d=3);
  * cross-checks automaton counts against a direct brute-force enumeration of
    memory-k walks in Z^3 for n <= 7, and against direct SAW counts;
  * replays the stored integer Collatz certificate q*A*w < p*w in exact
    integer arithmetic on the freshly rebuilt k=12 automaton;
  * re-derives the directed decimal floor of atanh(q/p) and compares it with
    the incumbent lower endpoint.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from collections import defaultdict, deque
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "bounds" / "finite_memory_saw.json"
DIM = 3
STEPS = tuple(
    tuple(s if i == a else 0 for i in range(DIM)) for a in range(DIM) for s in (1, -1)
)
SYMS = tuple(
    (p, s)
    for p in itertools.permutations(range(DIM))
    for s in itertools.product((1, -1), repeat=DIM)
)
INCUMBENT = "0.2122159753270231627267517174278577806020"
FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


# ------------------------------------------------------------ fresh automaton
def canon(path: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    origin = path[0]
    best = None
    for perm, signs in SYMS:
        image = tuple(
            tuple(signs[i] * (point[perm[i]] - origin[perm[i]]) for i in range(DIM))
            for point in path
        )
        if best is None or image < best:
            best = image
    assert best is not None
    return best


def can_close(path: tuple[tuple[int, ...], ...], memory: int) -> bool:
    """Fresh implementation: iterative depth-first search with L1 pruning."""

    budget = memory - (len(path) - 1)
    if budget < 1:
        return False
    start, end = path[0], path[-1]
    interior = set(path[1:-1])
    stack = [(end, 0)]
    seen = {(end, 0)}
    while stack:
        node, used = stack.pop()
        for step in STEPS:
            nxt = tuple(node[i] + step[i] for i in range(DIM))
            if nxt == start:
                return True
            spent = used + 1
            l1 = sum(abs(nxt[i] - start[i]) for i in range(DIM))
            if spent + l1 > budget or nxt in interior:
                continue
            key = (nxt, spent)
            if key not in seen:
                seen.add(key)
                stack.append(key)
    return False


CLOSE_CACHE: dict[tuple[tuple[tuple[int, ...], ...], int], bool] = {}


def shrink(path: tuple[tuple[int, ...], ...], memory: int) -> tuple[tuple[int, ...], ...]:
    while True:
        key = (path, memory)
        verdict = CLOSE_CACHE.get(key)
        if verdict is None:
            verdict = can_close(path, memory)
            CLOSE_CACHE[key] = verdict
        if verdict:
            return canon(path)
        path = path[1:]


def build(memory: int) -> tuple[list[tuple[tuple[int, ...], ...]], list[dict[int, int]]]:
    origin = (0,) * DIM
    seed = canon((origin, STEPS[0]))
    states = [seed]
    index = {seed: 0}
    rows: list[dict[int, int]] = []
    queue = deque([seed])
    while queue:
        state = queue.popleft()
        row: dict[int, int] = defaultdict(int)
        for step in STEPS:
            nxt = tuple(state[-1][i] + step[i] for i in range(DIM))
            if nxt in state:
                continue
            target = shrink(state + (nxt,), memory)
            slot = index.get(target)
            if slot is None:
                slot = len(states)
                index[target] = slot
                states.append(target)
                queue.append(target)
            row[slot] += 1
        rows.append(dict(row))
    return states, rows


def digest_states(states) -> str:
    import hashlib

    payload = ";".join(
        "/".join(",".join(str(c) for c in point) for point in state) for state in states
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def digest_rows(rows) -> str:
    import hashlib

    payload = ";".join(
        ",".join(f"{t}:{m}" for t, m in sorted(row.items())) for row in rows
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def automaton_counts(rows: list[dict[int, int]], length: int) -> list[int]:
    """Number of memory-k walks of length 1..length (6 first steps -> seed)."""

    weights = {0: 6}
    out = [6]
    for _ in range(length - 1):
        nxt: dict[int, int] = defaultdict(int)
        for state, count in weights.items():
            for target, multiplicity in rows[state].items():
                nxt[target] += count * multiplicity
        weights = dict(nxt)
        out.append(sum(weights.values()))
    return out


# --------------------------------------------------------- direct enumerations
def direct_memory_counts(memory: int, length: int) -> list[int]:
    """Brute-force count of walks with no revisit within `memory` steps."""

    counts = [0] * length
    origin = (0,) * DIM

    def recurse(path: list[tuple[int, ...]], depth: int) -> None:
        if depth == length:
            return
        window = path[max(0, len(path) - memory):]
        for step in STEPS:
            nxt = tuple(path[-1][i] + step[i] for i in range(DIM))
            if nxt in window:
                continue
            counts[depth] += 1
            path.append(nxt)
            recurse(path, depth + 1)
            path.pop()

    recurse([origin], 0)
    return counts


def direct_saw_counts(length: int) -> list[int]:
    counts = [0] * length
    origin = (0,) * DIM

    def recurse(path: set, last: tuple[int, ...], depth: int) -> None:
        if depth == length:
            return
        for step in STEPS:
            nxt = tuple(last[i] + step[i] for i in range(DIM))
            if nxt in path:
                continue
            counts[depth] += 1
            path.add(nxt)
            recurse(path, nxt, depth + 1)
            path.remove(nxt)

    recurse({origin}, origin, 0)
    return counts


def main() -> int:
    began = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    stored = {int(r["memory"]): r for r in artifact["data"]["memories"]}
    target = artifact["data"]["target"]
    target_memory = int(target["memory"])

    # 1. rebuild every automaton and compare structure digests
    rebuilt = {}
    for memory in sorted(stored):
        states, rows = build(memory)
        rebuilt[memory] = (states, rows)
        check(
            f"k={memory} rebuilt structure matches artifact",
            len(states) == int(stored[memory]["state_count"])
            and digest_states(states) == stored[memory]["state_sha256"]
            and digest_rows(rows) == stored[memory]["transition_sha256"],
            f"states={len(states)}",
        )

    # 2. exact k=4 characteristic polynomial versus the published closed form
    states4, rows4 = rebuilt[4]
    dense = [[rows4[i].get(j, 0) for j in range(len(states4))] for i in range(len(states4))]
    # charpoly of a 3x3 integer matrix, computed by hand
    a, b, c = dense[0], dense[1], dense[2]
    trace = a[0] + b[1] + c[2]
    minors = (
        b[1] * c[2] - b[2] * c[1]
        + a[0] * c[2] - a[2] * c[0]
        + a[0] * b[1] - a[1] * b[0]
    )
    det = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    check(
        "k=4 charpoly equals x^3-4x^2-4x-1 (Poenitz--Tittmann d=3 closed form)",
        len(states4) == 3 and trace == 4 and minors == -4 and det == 1,
        f"trace={trace}, second={minors}, det={det}",
    )

    # 3. automaton counts versus direct enumeration (n <= 7) and SAW domination
    saw = direct_saw_counts(7)
    known_saw = [6, 30, 150, 726, 3534, 16926, 81390]
    check("direct SAW counts match published c_1..c_7", saw == known_saw, str(saw))
    for memory in (4, 6):
        counts_auto = automaton_counts(rebuilt[memory][1], 7)
        counts_direct = direct_memory_counts(memory, 7)
        check(
            f"k={memory} automaton counts equal direct memory-walk counts, n<=7",
            counts_auto == counts_direct,
            str(counts_auto),
        )
        check(
            f"k={memory} counts dominate SAW counts",
            all(cs <= ca for cs, ca in zip(saw, counts_auto, strict=True)),
        )
    counts12 = automaton_counts(rebuilt[target_memory][1], 7)
    check(
        f"k={target_memory} counts equal SAW counts through n=7 (no loop that short)",
        counts12 == saw,
        str(counts12),
    )

    # 4. exact integer replay of the stored Collatz certificate on fresh rows
    weights = [int(w) for w in target["weights"]]
    numerator = int(target["exact_upper"]["numerator"])
    denominator = int(target["exact_upper"]["denominator"])
    states_t, rows_t = rebuilt[target_memory]
    check(
        "stored weight vector length matches rebuilt state count",
        len(weights) == len(states_t) and min(weights) >= 1,
    )
    residual_min = None
    strict = True
    for source, row in enumerate(rows_t):
        lhs = numerator * weights[source]
        rhs = denominator * sum(m * weights[t] for t, m in row.items())
        r = lhs - rhs
        if residual_min is None or r < residual_min:
            residual_min = r
        if r <= 0:
            strict = False
    check(
        f"integer Collatz replay: {numerator}*w > {denominator}*A*w componentwise",
        strict and residual_min == int(target["exact_upper"]["minimum_residual"]),
        f"minimum residual {residual_min}",
    )

    # 5. directed decimal floor of atanh(q/p) and incumbent comparison
    mp.mp.dps = 200
    value_hi = mp.atanh(mp.mpf(denominator) / numerator)
    mp.mp.dps = 120
    value_lo = mp.atanh(mp.mpf(denominator) / numerator)
    agree = mp.fabs(value_hi - value_lo) < mp.mpf(10) ** (-100)
    floor60 = mp.floor(value_lo * mp.mpf(10) ** 60) / mp.mpf(10) ** 60
    stored_floor = artifact["data"]["theorem"]["directed_floor"]
    check("atanh evaluation agrees at 120 and 200 dps", bool(agree))
    check(
        "stored directed floor matches recomputation",
        mp.nstr(floor60, 62) == stored_floor,
        stored_floor,
    )
    check(
        "new lower endpoint strictly exceeds the incumbent",
        floor60 > mp.mpf(INCUMBENT),
        f"{stored_floor} > {INCUMBENT}",
    )

    cpu = time.process_time() - began
    check("verifier process-time budget", cpu < 1200.0, f"{cpu:.1f}s")
    if FAILS:
        print(f"FAIL: {len(FAILS)} check(s) failed: {FAILS}")
        return 1
    print(f"PASS: independent finite-memory SAW verification complete ({cpu:.1f}s CPU)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
