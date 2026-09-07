"""Bounded source-admission probe, not a circuit-distance or fault-model solver.

Run against the registered campaign's immutable input snapshot. A nonzero
exit records an unadmitted source; it does not refute the d=11 conjecture.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    inputs = args.run / "inputs"
    prereg = json.loads((inputs / "preregistration.json").read_text())
    spec = json.loads((inputs / "figure_transcription.json").read_text())
    limits = prereg["limits"]
    resource.setrlimit(resource.RLIMIT_CPU, (limits["total_assay_wall_seconds"],) * 2)
    pins = prereg["source_pins"]
    source_files = {
        "2603.05481v1.tar.gz": pins["arxiv_source_sha256"],
        "STIM_LRCircuits.zip": pins["zenodo_archive_sha256"],
        "gross-code.pdf": pins["figure10_sha256"],
        "three-color-residuals.pdf": pins["figure11_sha256"],
        "three-color-circuits-compressed.pdf": pins["figure12_sha256"],
    }
    hash_checks = {
        name: hashlib.sha256((inputs / name).read_bytes()).hexdigest() == expected
        for name, expected in source_files.items()
    }
    if not all(hash_checks.values()):
        raise ValueError(f"Source hash mismatch: {hash_checks}")
    with zipfile.ZipFile(inputs / "STIM_LRCircuits.zip") as archive:
        inventory = [
            {"path": info.filename, "bytes": info.file_size,
             "sha256": hashlib.sha256(archive.read(info)).hexdigest()}
            for info in archive.infolist() if not info.is_dir()
        ]
    width, height = spec["width"], spec["height"]
    data = {(x, y) for y in range(height) for x in range(width) if (x + y) % 2}
    checks = []
    for style, geometry in spec["geometry"].items():
        offsets = {"N": (0, -1), "E": (1, 0), "S": (0, 1), "W": (-1, 0),
                   "NE": geometry["NE"], "SW": geometry["SW"]}
        for y in range(geometry["y_parity"], height, 2):
            for x in range(geometry["x_parity"], width, 2):
                colour = spec["colour_order"][(x // 2 + y // 2) % 3]
                edge_times = dict(zip(spec["edge_order"], spec["figure12_times"][style][colour], strict=True))
                gates = {
                    ((x + dx) % width, (y + dy) % height): edge_times[edge]
                    for edge, (dx, dy) in offsets.items()
                }
                checks.append({"style": style, "position": (x, y), "colour": colour, "gates": gates})
    filled = [check for check in checks if check["style"] == "filled"]
    hollow = [check for check in checks if check["style"] == "hollow"]
    odd_overlap_pairs = sum(
        len(a["gates"].keys() & b["gates"].keys()) % 2
        for a in filled for b in hollow
    )
    graph_checks = {
        "144_data": len(data) == 144,
        "72_checks_each": len(filled) == len(hollow) == 72,
        "all_weight_six": all(len(check["gates"]) == 6 for check in checks),
        "data_only_support": all(set(check["gates"]) <= data for check in checks),
        "CSS_commutation": odd_overlap_pairs == 0,
        "six_distinct_CNOT_times_per_ancilla": all(len(set(check["gates"].values())) == 6 for check in checks),
    }
    phases = []
    for phase in range(8):
        occupancy = defaultdict(list)
        for check in checks:
            offset = phase if check["style"] == "hollow" else 0
            for qubit, tick in check["gates"].items():
                occupancy[(qubit, (tick + offset) % 8)].append(
                    {"style": check["style"], "position": check["position"], "colour": check["colour"]})
        collisions = [(slot, owners) for slot, owners in occupancy.items() if len(owners) > 1]
        interleaved = []
        for a in filled:
            for b in hollow:
                overlap = a["gates"].keys() & b["gates"].keys()
                if not overlap:
                    continue
                for cycle in (-1, 0, 1):
                    diffs = [a["gates"][q] - (b["gates"][q] + phase + cycle * 8) for q in overlap]
                    if min(diffs) < 0 < max(diffs):
                        interleaved.append({"filled": a["position"], "hollow": b["position"], "relative_cycle": cycle})
        phases.append({
            "hollow_relative_origin": phase,
            "data_time_collision_slots": len(collisions),
            "first_collision": collisions[0] if collisions else None,
            "interleaved_check_pairs_with_cycle": len(interleaved),
            "first_interleaving": interleaved[0] if interleaved else None,
            "CNOT_only_feasible": not collisions and not interleaved and all(graph_checks.values()),
        })
    usage = resource.getrusage(resource.RUSAGE_SELF)
    report = {
        "status": "NOT_ADMITTED",
        "scope": "Literal Fig10/12 transcription and eight periodic relative-phase probes only; not the authors' validated circuit",
        "source_hash_checks": hash_checks,
        "archive_files": len(inventory),
        "archive_families": dict(Counter(Path(entry["path"]).name.split("_")[0] for entry in inventory)),
        "archive_gross_named_files": [entry["path"] for entry in inventory if any(token in entry["path"].lower() for token in ("gross", "bb144", "144,12,12"))],
        "graph_checks": graph_checks,
        "odd_XZ_overlap_pairs": odd_overlap_pairs,
        "periodic_CNOT_phase_probes": phases,
        "feasible_CNOT_only_phases": [p["hollow_relative_origin"] for p in phases if p["CNOT_only_feasible"]],
        "unresolved_source_obligations": spec["unresolved"],
        "fault_model_constructed": False,
        "direct_circuit_witness_replays": 0,
        "distance_solver_calls": 0,
        "distance_claim": None,
        "decision": "STOP distance search; resolve the exact timed-circuit/source contract first. No substitute circuit or boundary choice admitted.",
        "usage": {"wall_seconds": time.monotonic() - started, "cpu_seconds": usage.ru_utime + usage.ru_stime, "max_rss_native": usage.ru_maxrss},
        "effective_limits": {"CPU_seconds": list(resource.getrlimit(resource.RLIMIT_CPU)), "memory_enforcement": "See external RSS controller record; no hard address-space cap", "threads": 1},
    }
    (args.run / "archive_inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    (args.run / "source_admission.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
