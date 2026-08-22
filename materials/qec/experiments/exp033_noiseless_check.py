"""Stim noiseless-semantics check for the persisted EXP-033 depth-9 witnesses.

Proposition C5's schedules were verified with 0 detector firings in 4000
noiseless shots; EXP-033's witnesses were only structurally verified.  This
closes that gap: build the one-ancilla memory circuit from each persisted slot
map and require zero detector firings and zero observable flips noiselessly.

Writes an auxiliary artifact next to the canonical one; never touches it.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import importlib.util

import stim  # noqa: E402
from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    pbb_supports_and_orbits, pure_z_logical_basis)
from qec_research.circuits.scheduling import generator_supports, slots_to_layers  # noqa: E402
from qec_research.circuits.mixed_stabilizer import CircuitSpec, build_memory_circuit  # noqa: E402
from qec_research.artifacts import canonical_route  # noqa: E402

CANONICAL = ROOT / "results" / "processed" / "exp033_unrestricted_three.json"
OUT = ROOT / "results" / "processed" / "exp033_noiseless_check.json"
SHOTS = 4000
ROUNDS = 12
SEED = 20260813
EXPECTED_INSTANCES = {"pbb144-01", "pbb144-03", "pbb144-11"}
EXPECTED_DEPTH = 9


def canonical_json_sha256(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()

_SPEC = importlib.util.spec_from_file_location(
    "exp031_depth_criterion_for_noiseless", ROOT / "experiments" / "exp031_depth_criterion.py")
_E31 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _E31
_SPEC.loader.exec_module(_E31)


def write_result(
    payload: dict, *, clean: bool, full_coverage: bool, stamp: str | None = None
) -> Path:
    """Atomically route clean/full evidence; never overwrite it with a failure."""
    route = canonical_route(clean=clean, full_coverage=full_coverage)
    destination = OUT
    if route != "canonical":
        stamp = stamp or time.strftime("%Y%m%dT%H%M%S", time.gmtime())
        destination = (
            ROOT / "results" / route / f"exp033_noiseless_check-{stamp}.json"
        )
    payload = {
        **payload,
        "artifact_route": route,
        "artifact_path": str(destination.relative_to(ROOT)),
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2))
    tmp.replace(destination)
    return destination


def main() -> None:
    t0 = time.time()
    payload = json.loads(CANONICAL.read_text())
    canonical_input_sha256 = hashlib.sha256(CANONICAL.read_bytes()).hexdigest()
    instances = _E31._pbb_instances()
    rows = []
    for result in payload["results"]:
        index = int(result["exp031_instance"].removeprefix("pbb144-"))
        code, _, _ = pbb_supports_and_orbits(instances[index]["spec"])
        n = code.H.shape[1] // 2
        slot = {(check, qubit): depth for check, qubit, depth in result["slot"]}
        depth = max(slot.values()) + 1
        circuit, _meta = build_memory_circuit(CircuitSpec(
            H=code.H, observables=pure_z_logical_basis(code), rounds=ROUNDS,
            p=0.0, basis="Z", layers=slots_to_layers(slot, depth)))
        detector, observable = circuit.compile_detector_sampler(seed=SEED).sample(
            SHOTS, separate_observables=True)
        firings = int(detector.sum())
        flips = int(observable.sum())
        rows.append({
            "exp031_instance": result["exp031_instance"],
            "label": result["label"],
            "depth": depth,
            "noiseless_shots": SHOTS,
            "certified_depth": int(result["certified_depth"]),
            "w_max": int(result["w_max"]),
            "slot_sha256": canonical_json_sha256(result["slot"]),
            "rounds": ROUNDS,
            "noiseless_detector_firings": firings,
            "noiseless_observable_flips": flips,
            "passes": firings == 0 and flips == 0,
        })
        print(f"{result['exp031_instance']} ({result['label']}): depth={depth} "
              f"firings={firings} flips={flips}", flush=True)

    full_coverage = (
        len(rows) == len(EXPECTED_INSTANCES)
        and {row["exp031_instance"] for row in rows} == EXPECTED_INSTANCES
    )
    clean = full_coverage and all(
        row["passes"]
        and row["depth"] == row["certified_depth"] == row["w_max"] == EXPECTED_DEPTH
        for row in rows
    )
    verdict = "POSITIVE" if clean else "NEGATIVE"
    out_payload = {
        "experiment": "EXP-033-noiseless-check",
        "canonical_input": str(CANONICAL.relative_to(ROOT)),
        "canonical_input_sha256": canonical_input_sha256,
        "protocol": {"shots": SHOTS, "rounds": ROUNDS, "seed": SEED, "p": 0.0,
                     "basis": "Z"},
        "results": rows,
        "verdict": verdict,
        "wall_s": round(time.time() - t0, 1),
    }
    destination = write_result(
        out_payload, clean=clean, full_coverage=full_coverage
    )
    print(f"wrote {destination.relative_to(ROOT)}: {verdict}", flush=True)


if __name__ == "__main__":
    main()
