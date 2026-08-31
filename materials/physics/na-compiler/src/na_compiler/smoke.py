#!/usr/bin/env python3
"""End-to-end smoke: tiny instance, ILP solve vs QMAP, gap table row.

Runs the gate-A1 pipeline (QMAP compile -> parse -> certified optimum per
transition -> gap) on one 3-qubit QFT instance, exercises gate-B enumeration
with correctness invariants, and prints a bounded runtime-scaling table for
planning in-warehouse campaign budgets.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone

from .architecture import architecture_for
from .circuits import qft_instance, reg3_instance
from .config import TIME_ATOM_TRANSFER_US, move_batch_duration_us
from .conflict import Move, conflict_graph
from .optimize import cert_min_batches, cert_min_duration
from .qmap_adapter import QmapRunResult, compile_with_qmap, extract_placement_trajectory


def qmap_transport_cost_from_events(prog) -> dict:
    """QMAP output -> transport cost under the pinned duration model.

    One batch = one @+move event; its duration uses the max Euclidean atom
    travel inside that batch; each load and each store event pays
    TIME_ATOM_TRANSFER_US per atom. Positions update per batch in order.
    """
    positions = dict(prog.atom_locations)
    batches = 0
    batch_dists = []
    n_load_store = 0
    for ev in prog.events:
        if ev.kind == "load":
            n_load_store += len(ev.atoms)
        elif ev.kind == "store":
            n_load_store += len(ev.atoms)
        elif ev.kind == "move":
            batches += 1
            dmax = 0.0
            for atom, xy in ev.targets.items():
                assert atom in positions, f"moved atom {atom} has no prior position"
                dmax = max(dmax, math.dist(positions[atom], xy))
            batch_dists.append(dmax)
            for atom, xy in ev.targets.items():
                positions[atom] = xy
    total = (
        sum(move_batch_duration_us(d) for d in batch_dists)
        + TIME_ATOM_TRANSFER_US * n_load_store
    )
    return {
        "n_batches": batches,
        "batch_dists": batch_dists,
        "n_load_store": n_load_store,
        "total_us": total,
    }


def transition_moves(start: dict, target: dict) -> list[Move]:
    """Build the Move list for one placement transition."""
    moves = []
    for idx, atom in enumerate(sorted(start)):
        if start[atom] != target[atom]:
            moves.append(
                Move(
                    atom=idx,
                    start=start[atom],
                    target=target[atom],
                    start_site=None,
                    target_site=None,
                )
            )
    return moves


def run_gate_a(instance, verbose: bool = False) -> dict:
    """Gate A1: certified per-transition optimum vs QMAP transport cost."""
    arch = architecture_for(instance.n_qubits, max(len(l) for l in instance.pairs))
    run: QmapRunResult = compile_with_qmap(
        arch.to_quimap_json(), instance.pairs, instance.n_qubits
    )
    qmap_cost = qmap_transport_cost_from_events(run.program)
    traj = extract_placement_trajectory(
        run.program, [f"atom{i}" for i in range(instance.n_qubits)]
    )

    rows = []
    total_cert_duration = 0.0
    for t_idx in range(len(traj) - 1):
        moves = transition_moves(traj[t_idx], traj[t_idx + 1])
        if not moves:
            continue
        t0 = time.perf_counter()
        chi, witness, chi_log = cert_min_batches(moves)
        t_chi = time.perf_counter() - t0
        t0 = time.perf_counter()
        dur = cert_min_duration(moves, max_batches=len(moves))
        t_dur = time.perf_counter() - t0
        tot = dur["duration_opt"] + TIME_ATOM_TRANSFER_US * 2 * len(moves)
        total_cert_duration += tot
        rows.append(
            {
                "transition": t_idx,
                "n_moves": len(moves),
                "chi": chi,
                "chi_log": chi_log,
                "chi_time_s": t_chi,
                "duration_opt_us": dur["duration_opt"],
                "mip_log": dur["log"],
                "mip_time_s": t_dur,
                "cert_total_us": tot,
            }
        )
        if verbose:
            print(rows[-1])
    assert total_cert_duration > 0, "gate A instance must have moving transitions"
    gap = (qmap_cost["total_us"] - total_cert_duration) / total_cert_duration
    return {
        "instance": instance.name,
        "qmap_cost_us": qmap_cost["total_us"],
        "qmap_stats": {
            "n_batches": qmap_cost["n_batches"],
            "n_load_store": qmap_cost["n_load_store"],
        },
        "cert_duration_us": total_cert_duration,
        "gap_rel": gap,
        "qmap_wall_s": run.qmap_wall_s,
        "transitions": rows,
        "naviz": run.program.raw_events_log,
    }


def run_gate_b(
    n_atoms: int,
    grid_side: int,
    max_depth: int | None = None,
    max_states: int | None = None,
) -> dict:
    """Exhaustive BFS enumeration with correctness invariants."""
    from .exhaustive import enumerate_reachable, hardest_states

    t0 = time.perf_counter()
    dist, frontier_sizes, visited, aborted = enumerate_reachable(
        n_atoms, grid_side, max_depth, max_states
    )
    elapsed = time.perf_counter() - t0
    # invariants: states are injective maps into the fixed grid
    nsites = grid_side * grid_side
    from itertools import islice

    for s in islice(dist, 2000):
        assert len(set(s)) == n_atoms and max(s) < nsites, f"bad state {s}"
    hard = hardest_states(dist, 5)
    return {
        "n_atoms": n_atoms,
        "grid_side": grid_side,
        "reachable_states": visited,
        "max_depth": max(d for _, d in dist.items()),
        "frontier_sizes": frontier_sizes,
        "elapsed_s": elapsed,
        "hardest": [{"state": list(s), "min_batches": d} for s, d in hard],
        "complete_enumeration": (not aborted),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    started = datetime.now(timezone.utc).isoformat()

    print(f"# na-compiler smoke run {started}")
    print(f"# python {platform.python_version()} on {platform.machine()}")

    inst = qft_instance(3)
    print(f"\n## gate A: {inst.name}, layers={[len(l) for l in inst.pairs]}")
    a1 = run_gate_a(inst, verbose=args.verbose)
    print(
        json.dumps(
            {k: v for k, v in a1.items() if k != "transitions"},
            indent=2,
            default=str,
        )
    )
    for r in a1["transitions"]:
        print(
            f"  transition {r['transition']}: moves={r['n_moves']} chi={r['chi']} "
            f"cert={r['cert_total_us']:.2f}us chi_log={r['chi_log'][:70]!r} "
            f"mip_log={r['mip_log']!r}"
        )

    print("\n## gate B:")
    for n, side in [(2, 2), (3, 2), (3, 3)]:
        b = run_gate_b(n, side)
        print(
            f"  n={n} grid={side}x{side}: states={b['reachable_states']} "
            f"maxdepth={b['max_depth']} time={b['elapsed_s']:.2f}s "
            f"complete={b['complete_enumeration']} hardest={b['hardest'][:2]}"
        )

    print("\n## scaling (gate B exhaustive BFS, bounded):")
    for n, side in [(3, 3), (4, 3), (5, 3)]:
        b = run_gate_b(n, side)
        print(
            f"  n={n} grid={side}x{side}: states={b['reachable_states']:<8} "
            f"time={b['elapsed_s']:>7.2f}s maxdepth={b['max_depth']} "
            f"complete={b['complete_enumeration']}"
        )
    # 4x4 probes bound: cap at 30k states so the full smoke remains wall-
    # bounded; cap-reached rows are labeled hardest-found (complete=False).
    for n in (4, 5, 6, 7):
        b = run_gate_b(n, 4, max_states=30_000)
        print(
            f"  n={n} grid=4x4: states={b['reachable_states']:<8} "
            f"time={b['elapsed_s']:>7.2f}s maxdepth={b['max_depth']} "
            f"complete={b['complete_enumeration']} (cap 30k)"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
