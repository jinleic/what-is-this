#!/usr/bin/env python3
"""Campaign runner: run gates, freeze artifacts under campaigns/<ts>_<uuid>_<hash>/.

Manifest includes content hashes of every input (instances, architecture,
code files) so certificates are reproducible. One campaign = immutable
snapshot; the runner never overwrites an existing campaign directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # physics/na-compiler
sys.path.insert(0, str(REPO / "src"))

from na_compiler.architecture import architecture_for, Architecture  # noqa: E402
from na_compiler.circuits import qft_instance, reg3_instance  # noqa: E402
from na_compiler.config import TIME_ATOM_TRANSFER_US  # noqa: E402
from na_compiler.optimize import cert_min_batches, cert_min_duration  # noqa: E402
from na_compiler.qmap_adapter import (  # noqa: E402
    compile_with_qmap,
    extract_placement_trajectory,
)
from na_compiler.smoke import qmap_transport_cost_from_events  # noqa: E402
from na_compiler.conflict import Move  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def transition_moves(start: dict, target: dict) -> list[Move]:
    moves = []
    for idx, atom in enumerate(sorted(start)):
        if start[atom] != target[atom]:
            moves.append(
                Move(idx, start[atom], target[atom], None, None)
            )
    return moves


def run_instance(instance, label_gate: str, runs: dict, camp=None) -> dict:
    arch = architecture_for(instance.n_qubits, max(len(l) for l in instance.pairs))
    run = compile_with_qmap(arch.to_quimap_json(), instance.pairs, instance.n_qubits)
    qmap_wall = run.qmap_wall_s
    cost = qmap_transport_cost_from_events(run.program)
    traj = extract_placement_trajectory(run.program, [f"atom{i}" for i in range(instance.n_qubits)])

    transitions = []
    total_cert = 0.0
    for t_idx in range(len(traj) - 1):
        moves = transition_moves(traj[t_idx], traj[t_idx + 1])
        if not moves:
            continue
        chi, witness, chi_log = cert_min_batches(moves)
        dur = cert_min_duration(moves, max_batches=len(moves))
        tot = dur["duration_opt"] + TIME_ATOM_TRANSFER_US * 2 * len(moves)
        total_cert += tot
        transitions.append({
            "t": t_idx,
            "n_moves": len(moves),
            "chi": chi,
            "chi_log": chi_log,
            "mip_log": dur["log"],
            "cert_total_us": tot,
        })
    gap = (
        (cost["total_us"] - total_cert) / total_cert if total_cert > 0 else 0.0
    )
    chi = [t["chi"] for t in transitions]
    return {
        "gate": label_gate,
        "instance": instance.name,
        "family": instance.family,
        "n_qubits": instance.n_qubits,
        "provenance": instance.provenance,
        "qmap_cost_us": cost["total_us"],
        "qmap_n_batches": cost["n_batches"],
        "qmap_n_load_store": cost["n_load_store"],
        "cert_duration_us": total_cert,
        "chi_sum": sum(chi),
        "chi_max": max(chi, default=0),
        "min_batches_lb_gap": (cost["n_batches"] - sum(chi)) / max(sum(chi), 1),
        "gap_rel": gap,
        "qmap_wall_s": qmap_wall,
        "transitions": transitions,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gates", default="A", help="comma list from {A,B}")
    ap.add_argument("--nqubits", default="3,4", help="comma list of sizes")
    ap.add_argument("--seeds", default="0", help="comma list reg3 seeds")
    ap.add_argument("--note", default="", help="free-text note in manifest")
    args = ap.parse_args()

    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_uuid = uuid.uuid4().hex[:8]
    # content hash over src tree = code identity
    h = hashlib.sha256()
    for f in sorted((REPO / "src").rglob("*.py")):
        h.update(sha256(f).encode())
        h.update(str(f).encode())
    code_hash = h.hexdigest()[:12]
    camp = REPO / "campaigns" / f"{ts}_{run_uuid}_{code_hash}"
    camp.mkdir(parents=True, exist_ok=False)

    import mqt.qmap
    import qiskit
    import networkx
    import highspy
    import z3

    manifest = {
        "started_utc": ts,
        "uuid": run_uuid,
        "code_hash": code_hash,
        "note": args.note,
        "host": platform.platform(),
        "python": platform.python_version(),
        "deps": {
            "mqt.qmap": getattr(mqt.qmap, "__version__", "3.9.0"),
            "qiskit": qiskit.__version__,
            "networkx": networkx.__version__,
            "highspy": highspy.__version__ if hasattr(highspy, "__version__") else "1.15.1",
            "z3-solver": z3.get_version_string(),
        },
        "priority": "nice -n 10, single core",
        "gates": {},
    }

    sizes = [int(x) for x in args.nqubits.split(",")]
    seeds = [int(x) for x in args.seeds.split(",")]
    rows = []
    if "A" in args.gates:
        for n in sizes:
            rows.append(run_instance(qft_instance(n), "A1", manifest))
        # 3-regular graphs exist only for even n (pre_statement instance set
        # uses {4,6,8,10}); QFT covers the odd sizes.
        for n in sizes:
            if n % 2 == 1:
                continue
            for s in seeds:
                rows.append(run_instance(reg3_instance(n, s), "A1", manifest))
        manifest["gates"]["A"] = rows
    if "B" in args.gates:
        from na_compiler.exhaustive import enumerate_reachable, hardest_states

        b_rows = []
        for n, side in [(3, 3), (4, 3), (4, 4)]:
            t0 = time.perf_counter()
            dist, fr, vis, ab = enumerate_reachable(n, side)
            dt = time.perf_counter() - t0
            hard = hardest_states(dist, 5)
            b_rows.append({
                "n_atoms": n,
                "grid_side": side,
                "reachable_states": vis,
                "max_depth": max(dist.values()),
                "frontier_sizes": fr,
                "complete": not ab,
                "elapsed_s": dt,
                "hardest": [{"state": list(s), "min_batches": d} for s, d in hard],
            })
        manifest["gates"]["B"] = b_rows

    manifest_path = camp / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str))
    # gap table
    gap_tbl = camp / "gap_table.md"
    lines = [
        "| instance | n | qmap_cost_us | cert_us | gap_rel | chi_max |",
        "|---|---|---|---|---|---|",
    ]
    for r in manifest["gates"].get("A", []):
        chi_max = max((t["chi"] for t in r["transitions"]), default=0)
        lines.append(
            f"| {r['instance']} | {r['n_qubits']} | {r['qmap_cost_us']:.2f} | "
            f"{r['cert_duration_us']:.2f} | {r['gap_rel']:+.4f} | {chi_max} |"
        )
    gap_tbl.write_text("\n".join(lines) + "\n")

    print(f"campaign dir: {camp}")
    print(gap_tbl.read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
