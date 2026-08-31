"""Gate-B COMPANION at elevated p (explicitly NOT the pinned gate).

The pinned Gate B compares beam8 with bp30+osd at p=1e-3, where both arms
produced 0 failures in 1e4 shots (UNDERDETERMINED). Resolving it at the
pinned p needs >=1e6 shots/arm. This companion instead measures the SAME
ratio at elevated p, where failures are resolvable in ~1e3 shots, to see
whether the equal-accuracy claim's direction is even plausible before
committing ~71 core-h.

Scope discipline:
  * This is a COMPANION measurement. It cannot pass or fail the pinned
    Gate B (different p), and the campaign is labelled accordingly.
  * The elevated-p circuits do not exist upstream. They are DERIVED from
    the committed IonQ p=1e-3 circuit by uniform rescale of every noise
    argument -- the same documented technique already used in-repo for
    the p=3e-4 point -- and every derived file is sha256-logged.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import stim

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qldpc_dec.bootstrap import ratio_ci, wilson_interval  # noqa: E402
from qldpc_dec.circuits import CIRCUIT_DIR, load_circuit  # noqa: E402
from qldpc_dec.dem_matrices import dem_to_matrices  # noqa: E402
from qldpc_dec.runner import Campaign  # noqa: E402
from qldpc_dec.run_gate import LER_NORMALIZATION_ROUNDS  # noqa: E402
from qldpc_dec.seeds import sampling_seed  # noqa: E402

BASE_P = 1e-3


def derived_circuit(p_target: float) -> tuple[stim.Circuit, str, Path]:
    """Uniformly rescale every noise argument of the committed p=1e-3 Z circuit."""
    base = load_circuit(BASE_P, "Z")
    text = str(base)
    factor = p_target / BASE_P

    def scale(m: re.Match) -> str:
        return f"({float(m.group(1)) * factor:.12g})"

    scaled = re.sub(r"\(([0-9.eE+-]+)\)", scale, text)
    circ = stim.Circuit(scaled)
    path = CIRCUIT_DIR / (
        f"BB_144_144_12_memory_Z_p{p_target:g}_sr12_derived_p1e-3_rescale.stim"
    )
    path.write_text(str(circ))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    return circ, sha, path


def arm_run(arm: str, circ: stim.Circuit, dem, shots: int, seed: int) -> dict:
    det, obs = circ.compile_detector_sampler(seed=seed).sample(
        shots, separate_observables=True
    )
    t0 = time.perf_counter()
    if arm == "bp30+osd":
        from qldpc_dec.bp_osd import bp_osd_harness

        preds = bp_osd_harness().decode(dem, det)
    elif arm == "beam8":
        from qldpc_dec.beam_search import BeamSearchBatchDecoder

        dec = BeamSearchBatchDecoder(
            dem, beam_width=8, initial_iters=30, iters_per_round=20,
            max_rounds=10, num_results=1,
        )
        preds = dec.decode_batch(det)
    else:
        raise ValueError(arm)
    wall = time.perf_counter() - t0
    fails = int((preds != obs).any(axis=1).sum())
    raw_rate = fails / shots
    raw_ci = wilson_interval(fails, shots)
    return {
        "arm": arm,
        "shots": shots,
        "failures": fails,
        "raw_logical_failure_rate": raw_rate,
        "raw_ler_ci95": list(raw_ci),
        "per_round_ler": raw_rate / LER_NORMALIZATION_ROUNDS,
        "per_round_ler_ci95": [
            endpoint / LER_NORMALIZATION_ROUNDS for endpoint in raw_ci
        ],
        "ler_normalization_rounds": LER_NORMALIZATION_ROUNDS,
        "wall_s": wall,
        "ms_per_shot": 1e3 * wall / shots,
        "seed": seed,
    }


def main(p_target: float, shots: int) -> None:
    camp = Campaign({"run": "gateB_companion_highp", "p": p_target, "shots": shots,
                     "base_seed": 20260830})
    camp.write_manifest()
    circ, sha, path = derived_circuit(p_target)
    dem = circ.detector_error_model(
        decompose_errors=True, ignore_decomposition_failures=True
    )
    results = {}
    for arm in ("bp30+osd", "beam8"):
        seed = sampling_seed(20260830, p_target, "Z", arm)
        r = arm_run(arm, circ, dem, shots, seed)
        camp.append_result(r)
        results[arm] = r
        print(f"{arm}: raw={r['raw_logical_failure_rate']:.5g}, "
              f"per-round LER={r['per_round_ler']:.5g} "
              f"({r['failures']}/{shots}) {r['ms_per_shot']:.0f} ms/shot")

    a, b = results["beam8"], results["bp30+osd"]
    lo, hi = ratio_ci(a["failures"], a["shots"], b["failures"], b["shots"])
    verdict = {
        "companion_p": p_target,
        "beam8_over_bposd_ler_ratio": (
            a["per_round_ler"] / b["per_round_ler"]
        ) if b["per_round_ler"] else None,
        "ratio_ci95": [lo, hi],
        "pinned_gate_band": [0.87, 1.15],
        "band_consistent_at_this_p": bool(lo <= 1.15 and hi >= 0.87),
        "scope": "COMPANION at elevated p; cannot pass/fail the pinned p=1e-3 Gate B",
        "derived_circuit_sha256": sha,
        "derived_circuit": path.name,
    }
    camp.close({"gate": "B-companion-highp", **verdict,
                "bposd": b, "beam8": a})
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 3e-3,
         int(sys.argv[2]) if len(sys.argv) > 2 else 2000)
