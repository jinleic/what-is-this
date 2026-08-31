"""Owner verification of NaCompilerExt's gate-B optimal-COST rows.

Independent of the agent's code path: re-implements witness replay and
cost recomputation directly from the FROZEN primitives
(`na_compiler.conflict.compatible` for legality, `config.move_batch_duration_us`
+ `TIME_ATOM_TRANSFER_US` for duration) and checks:

  1. each witness batch is pairwise-legal and moves to free, distinct sites;
  2. applying the witness reproduces the claimed intermediate states and
     lands exactly on the claimed target permutation;
  3. the witness cost under the pinned model equals the claimed opt_cost_us;
  4. the batch count matches the owner-verified BFS diameter for that n
     (4 for n=5, 5 for n=6) — i.e. the witness is batch-optimal.

Cost-OPTIMALITY of the value itself (that no cheaper 4/5-batch schedule
exists) is the agent's A* claim and is labelled [DERIVED] unless the
exact recheck below (small n only) is run.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from na_compiler.conflict import Move, compatible  # noqa: E402
from na_compiler.config import (  # noqa: E402
    STORAGE_SITE_SEPARATION,
    TIME_ATOM_TRANSFER_US,
    move_batch_duration_us,
)

GRID = 4
SEP_X, SEP_Y = STORAGE_SITE_SEPARATION


def site_xy(site: int) -> tuple[float, float]:
    return ((site % GRID) * SEP_X, (site // GRID) * SEP_Y)


def batch_legal(batch: list[list[int]], state: tuple[int, ...]) -> tuple[bool, str]:
    """Pairwise legality + free distinct targets, against the frozen relation.

    Pair semantics (as used by the producing agent and confirmed by exact
    witness replay): (atom_index, target_site).
    """
    moves = []
    for atom, dst in batch:
        if not (0 <= atom < len(state)):
            return False, f"atom index {atom} out of range for {state}"
        src = state[atom]
        moves.append(Move(atom, site_xy(src), site_xy(dst), src, dst))
    dsts = [d for _a, d in batch]
    if len(set(dsts)) != len(dsts):
        return False, "duplicate targets in batch"
    movers = {a for a, _d in batch}
    stationary = {s for i, s in enumerate(state) if i not in movers}
    for d in dsts:
        if d in stationary:
            return False, f"target {d} occupied by a non-moving atom"
    for i in range(len(moves)):
        for j in range(i + 1, len(moves)):
            if not compatible(moves[i], moves[j]):
                return False, f"moves {batch[i]} and {batch[j]} incompatible"
    return True, "ok"


def apply_batch(state: tuple[int, ...], batch: list[list[int]]) -> tuple[int, ...]:
    s = list(state)
    for atom, dst in batch:
        s[atom] = dst
    return tuple(s)


def batch_cost_us(batch: list[list[int]], state: tuple[int, ...]) -> float:
    dmax = 0.0
    for atom, dst in batch:
        x0, y0 = site_xy(state[atom])
        x1, y1 = site_xy(dst)
        d = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5
        dmax = max(dmax, d)
    return move_batch_duration_us(dmax) + 2.0 * TIME_ATOM_TRANSFER_US * len(batch)



def verify(instance: dict, expected_batches: int) -> dict:
    states = [tuple(s) for s in instance["witness_states"]]
    batches = instance["witness_batches"]
    out = {"instance": instance["instance"], "checks": {}}

    st = states[0]
    legal_all, replay_ok, total = True, True, 0.0
    for i, b in enumerate(batches):
        ok, why = batch_legal(b, st)
        if not ok:
            legal_all = False
            out["checks"][f"batch{i}_legality"] = why
        total += batch_cost_us(b, st)
        st = apply_batch(st, b)
        if st != states[i + 1]:
            replay_ok = False
            out["checks"][f"batch{i}_state"] = f"got {st}, witness says {states[i+1]}"

    target = tuple(
        int(x) for x in instance["instance"].rsplit("_", 1)[-1].split("-")
    )
    out["checks"]["all_batches_legal"] = legal_all
    out["checks"]["replay_matches_witness_states"] = replay_ok
    out["checks"]["lands_on_target"] = st == target
    out["final_state"] = list(st)
    out["target_state"] = list(target)
    out["recomputed_cost_us"] = total
    out["claimed_cost_us"] = instance["opt_cost_us"]
    out["cost_abs_diff_us"] = abs(total - instance["opt_cost_us"])
    out["cost_matches_1e-9"] = out["cost_abs_diff_us"] < 1e-9
    out["batches"] = len(batches)
    out["batch_count_matches_bfs_diameter"] = len(batches) == expected_batches
    return out


if __name__ == "__main__":
    payload = json.load(open(sys.argv[1]))
    expected = {5: 4, 6: 5}
    results = []
    for inst in payload["instances_closed"]:
        n = len(inst["witness_states"][0])
        results.append(verify(inst, expected[n]))
    print(json.dumps(results, indent=1))
