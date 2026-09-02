"""End-to-end gate smoke driver + campaign entry points.

Usage:
    python -m qldpc_dec.run_gate smoke             # 1e3-shot end-to-end proof
    python -m qldpc_dec.run_gate gateA_accuracy    # BP+OSD LERs at 3 pts
    python -m qldpc_dec.run_gate gateA_timing      # per-shot timing, Table II/III arm
    python -m qldpc_dec.run_gate gateB             # beam8 LER ratio at p=1e-3
    python -m qldpc_dec.run_gate gateC_floor       # ensemble-NMS vs 1.93e-9 floor

All runs: nice-bounded process, fixed seeds (seeds.py), snapshots into
campaigns/<UTC>_<uuid8>_<hash12>/ via runner.Campaign. The smoke run is a
REHEARSAL and stores under campaigns-smoke/ per physics/README.md.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import shutil
import time
from pathlib import Path

import numpy as np
import stim

from .bootstrap import bootstrap_ci_binomial, ratio_ci, wilson_interval
from .circuits import circuit_sha, dem_sha, load_circuit, load_dem
from .dem_matrices import dem_to_matrices
from .runner import Campaign
from .seeds import sampling_seed

# Official IonQ `simulation_functions.py`: per-basis LER is the any-logical
# shot-failure probability divided by the code distance / syndrome rounds.
LER_NORMALIZATION_ROUNDS = 12
BPOSD_PRIOR_MODEL = "merged DEM heterogeneous error_channel"
BPOSD_REFERENCE = {
    "package": "stimbposd==0.1.0",
    "commit": "7921f5eb1b358ff616f9822280c9961e83df06cb",
}
_SOURCE_DIR = Path(__file__).resolve().parent
_BPOSD_SNAPSHOT_SOURCES = (
    "run_gate.py",
    "bp_osd.py",
    "dem_matrices.py",
    "bootstrap.py",
    "circuits.py",
    "seeds.py",
    "runner.py",
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _write_bposd_manifest(campaign: Campaign, extra: dict | None = None) -> dict:
    source_hashes: dict[str, str] = {}
    for name in _BPOSD_SNAPSHOT_SOURCES:
        source = _SOURCE_DIR / name
        snapshot = campaign.dir / f"source_qldpc_dec_{name}"
        shutil.copy2(source, snapshot)
        source_hashes[snapshot.name] = _sha256_file(snapshot)
    pre_statement = _SOURCE_DIR.parent.parent / "pre_statement.md"
    pre_statement_snapshot = campaign.dir / "pre_statement.md"
    shutil.copy2(pre_statement, pre_statement_snapshot)
    manifest_extra = {
        "source_sha256": source_hashes,
        "pre_statement_sha256": _sha256_file(pre_statement_snapshot),
        "bposd_prior_model": BPOSD_PRIOR_MODEL,
        "bposd_reference": BPOSD_REFERENCE,
        "versions": {
            "numpy": np.__version__,
            "stim": stim.__version__,
            "ldpc": importlib.metadata.version("ldpc"),
        },
    }
    if extra:
        manifest_extra.update(extra)
    return campaign.write_manifest(manifest_extra)


def _within_factor_two_ratio(ratio: float) -> bool:
    return bool(np.isfinite(ratio) and 0.5 <= ratio <= 2.0)


def _bposd_prior_check(
    dem: stim.DetectorErrorModel, decoder: object
) -> dict[str, object]:
    _H, _A, priors = dem_to_matrices(dem, merge=True)
    priors_f64 = np.asarray(priors, dtype="<f8")
    configured = np.asarray(getattr(decoder, "channel_probs"), dtype="<f8")
    configured_exact = bool(np.array_equal(configured, priors_f64))
    if not configured_exact:
        raise AssertionError("BP+OSD configured channel differs from merged DEM priors")
    return {
        "columns": int(priors_f64.size),
        "distinct_probabilities": int(np.unique(priors_f64).size),
        "minimum": float(priors_f64.min()),
        "maximum": float(priors_f64.max()),
        "sha256_float64": hashlib.sha256(priors_f64.tobytes()).hexdigest(),
        "configured_channel_exact": configured_exact,
    }


def sample_shots(circuit: stim.Circuit, shots: int, seed: int):
    sampler = circuit.compile_detector_sampler(seed=seed)
    det, obs = sampler.sample(shots, separate_observables=True)
    return det, obs


def failure_stats(predictions: np.ndarray, observables: np.ndarray) -> tuple[int, int]:
    """Any-observable-mismatch = failure; returns (failures, shots)."""
    mism = np.any(predictions != observables, axis=1)
    return int(mism.sum()), int(len(mism))


def run_point(
    arm: str,
    p: float,
    basis: str,
    shots: int,
    base_seed: int,
    campaign: Campaign,
    decoder_kwargs: dict | None = None,
    time_it: bool = False,
) -> dict:
    """Decode one point; ``ler`` follows IonQ's per-round raw/12 convention."""
    decoder_kwargs = decoder_kwargs or {}
    circuit = load_circuit(p, basis)
    dem = load_dem(p, basis)
    seed = sampling_seed(base_seed, p, basis, arm)
    det, obs = sample_shots(circuit, shots, seed)

    t0 = time.perf_counter()
    if arm == "bp30+osd":
        from .bp_osd import bp_osd_harness

        harness = bp_osd_harness(**decoder_kwargs)
        preds = harness.decode(dem, det)
        times = None
    elif arm == "bp30+osd10k":
        from .bp_osd import bp_osd_harness

        harness = bp_osd_harness(max_iter=10000, **decoder_kwargs)
        preds = harness.decode(dem, det)
        times = None
    elif arm == "ensemble-nms24":
        from .nms import NmsEnsembleDecoder

        nms_kwargs = dict(ensemble_size=24)
        nms_kwargs.update(decoder_kwargs)
        dec = NmsEnsembleDecoder(dem, **nms_kwargs)
        preds = dec.decode_batch(det)
        times = None
        comp_stats = dict(dec.stats)
    elif arm.startswith("beam"):
        from .beam_search import BeamSearchBatchDecoder

        cfg = {
            "beam8": dict(beam_width=8, initial_iters=30, iters_per_round=20, max_rounds=10, num_results=1),
            "beam32": dict(beam_width=32, initial_iters=40, iters_per_round=30, max_rounds=10, num_results=1),
            "beam64": dict(beam_width=64, initial_iters=40, iters_per_round=30, max_rounds=20, num_results=1),
        }[arm]
        cfg.update(decoder_kwargs)
        dec = BeamSearchBatchDecoder(dem, **cfg)
        preds = dec.decode_batch(det)
        times = None
    else:
        raise ValueError(f"unknown arm {arm}")
    t1 = time.perf_counter()

    fails, n = failure_stats(preds, obs)
    raw_ler = fails / n
    raw_ci = wilson_interval(fails, n)
    result = {
        "arm": arm,
        "p": p,
        "basis": basis,
        "shots": n,
        "failures": fails,
        "raw_logical_failure_rate": raw_ler,
        "raw_ler_ci95": list(raw_ci),
        "ler": raw_ler / LER_NORMALIZATION_ROUNDS,
        "ler_ci95": [x / LER_NORMALIZATION_ROUNDS for x in raw_ci],
        "ler_normalization_rounds": LER_NORMALIZATION_ROUNDS,
        "sampling_seed": seed,
        "wall_s": t1 - t0,
        "circuit_sha256": circuit_sha(p, basis),
        "dem_sha256": dem_sha(p, basis),
    }
    if arm.startswith("bp30+osd"):
        result["bposd_prior_model"] = BPOSD_PRIOR_MODEL
        result["bposd_reference"] = BPOSD_REFERENCE
    if arm == "ensemble-nms24":
        result["ensemble_stats"] = comp_stats
    campaign.append_result(result)
    return result


def run_timing_point(
    p: float,
    basis: str,
    num_timed: int,
    base_seed: int,
    campaign: Campaign,
) -> dict:
    """Gate A timing arm: IonQ decoding_time protocol.

    10k single-decode() timings of ONE compiled decoder instance on the
    SAME pre-sampled shots (perf_counter each), reported as mean and
    p99.9 (their tail_index = int(0.999*n)-1).
    """
    circuit = load_circuit(p, basis)
    dem = load_dem(p, basis)
    seed = sampling_seed(base_seed, p, basis, "timing")
    det, _ = sample_shots(circuit, num_timed, seed)

    from .bp_osd import make_bp_osd_decoder

    dec = make_bp_osd_decoder(dem, max_iter=30)  # ms / osd_cs / order 10
    prior_check = _bposd_prior_check(dem, dec)
    times = np.empty(num_timed)
    for i in range(num_timed):
        t0 = time.perf_counter()
        dec.decode(det[i].astype(np.uint8))
        times[i] = time.perf_counter() - t0
    times_sorted = np.sort(times)
    tail_index = int(0.999 * num_timed) - 1
    result = {
        "arm": "bp30+osd-timing",
        "p": p,
        "basis": basis,
        "num_timed": num_timed,
        "mean_ms": float(times.mean() * 1e3),
        "p999_ms": float(times_sorted[tail_index] * 1e3),
        "p50_ms": float(np.median(times) * 1e3),
        "p99_ms": float(times_sorted[int(0.99 * num_timed) - 1] * 1e3),
        "sampling_seed": seed,
        "protocol": "single-decode perf_counter, one compiled instance (IonQ decoding_time)",
        "circuit_sha256": circuit_sha(p, basis),
        "dem_sha256": dem_sha(p, basis),
        "bposd_prior_model": BPOSD_PRIOR_MODEL,
        "bposd_reference": BPOSD_REFERENCE,
        "bposd_prior_check": prior_check,
    }
    campaign.append_result(result)
    return result


SMOKE_SHOTS = 1000


def cmd_smoke(base_seed: int = 20260829):
    """1e3-shot end-to-end rehearsal under campaigns-smoke/ (not a campaign)."""
    smoke_root = Path(__file__).resolve().parent.parent.parent / "campaigns-smoke"
    smoke_root.mkdir(exist_ok=True)
    camp = Campaign({"run": "smoke", "shots": SMOKE_SHOTS}, root=smoke_root)
    _write_bposd_manifest(
        camp,
        {"rehearsal": True, "note": "1e3-shot end-to-end proof, NOT a campaign"},
    )
    print("smoke campaign dir:", camp.dir)

    t_all0 = time.perf_counter()
    # gate A accuracy path at one point
    r1 = run_point("bp30+osd", 1e-3, "Z", SMOKE_SHOTS, base_seed, camp)
    print(f"[bp30+osd  p=1e-3 Z] LER={r1['ler']:.4g} CI95={r1['ler_ci95']} wall={r1['wall_s']:.1f}s "
          f"({r1['wall_s']/SMOKE_SHOTS*1e3:.2f} ms/shot)")

    # gate B path at one point
    r2 = run_point("beam8", 1e-3, "Z", SMOKE_SHOTS, base_seed, camp)
    print(f"[beam8     p=1e-3 Z] LER={r2['ler']:.4g} CI95={r2['ler_ci95']} wall={r2['wall_s']:.1f}s "
          f"({r2['wall_s']/SMOKE_SHOTS*1e3:.2f} ms/shot)")

    # gate C path at one point (tiny ensemble budget for the rehearsal)
    r3 = run_point("ensemble-nms24", 1e-3, "Z", 25, base_seed, camp,
                   decoder_kwargs={"ensemble_size": 2, "max_iters": 120})
    print(f"[nms ens2  p=1e-3 Z] 25 shots wall={r3['wall_s']:.1f}s "
          f"({r3['wall_s']/25*1e3:.1f} ms/shot) stats={r3.get('ensemble_stats')}")

    # timing path with a small timed count
    rt = run_timing_point(1e-3, "Z", 50, base_seed, camp)
    print(f"[bp30+osd timing x50] mean={rt['mean_ms']:.2f}ms p999={rt['p999_ms']:.2f}ms")

    t_all1 = time.perf_counter()
    summary = {
        "kind": "smoke-rehearsal",
        "total_wall_s": t_all1 - t_all0,
        "ms_per_shot_bp30osd": r1["wall_s"] / SMOKE_SHOTS * 1e3,
        "ms_per_shot_beam8": r2["wall_s"] / SMOKE_SHOTS * 1e3,
        "ms_per_shot_nmsens2": r3["wall_s"] / 25 * 1e3,
        "bp30osd_timing": {k: rt[k] for k in ("mean_ms", "p999_ms", "p50_ms", "p99_ms")},
        "gateA_cost_estimate_full": {
            "shots_per_point": 100_000,
            "points": 6,
            "single_core_hours": (r1["wall_s"] / SMOKE_SHOTS) * 100_000 * 6 / 3600,
        },
        "falsification_note": "smoke proves pipeline-end-to-end only; no gate verdicts",
    }
    camp.close(summary)
    print(json.dumps(summary, indent=2))


def cmd_gate_a_accuracy(shots: int, base_seed: int = 20260829):
    camp = Campaign({
        "run": "gateA_accuracy",
        "shots": shots,
        "base_seed": base_seed,
        "bposd_prior_model": BPOSD_PRIOR_MODEL,
    })
    _write_bposd_manifest(camp)
    results = []
    for p in (3e-4, 5e-4, 1e-3):
        for basis in ("X", "Z"):
            r = run_point("bp30+osd", p, basis, shots, base_seed, camp)
            print(f"p={p:g} {basis}: LER={r['ler']:.4g} CI={r['ler_ci95']}")
            results.append(r)
    xz = {}
    for p in (3e-4, 5e-4, 1e-3):
        lx = next(r for r in results if r["p"] == p and r["basis"] == "X")
        lz = next(r for r in results if r["p"] == p and r["basis"] == "Z")
        xz[f"{p:g}"] = {"total_ler": lx["ler"] + lz["ler"]}
    camp.close({"gate": "A-accuracy-companion", "total_ler_by_p": xz})


def cmd_gate_a_timing(num_timed: int, base_seed: int = 20260829):
    camp = Campaign({
        "run": "gateA_timing",
        "num_timed": num_timed,
        "base_seed": base_seed,
        "bposd_prior_model": BPOSD_PRIOR_MODEL,
    })
    _write_bposd_manifest(camp)
    out = {}
    for p in (1e-3, 5e-4):
        r = run_timing_point(p, "Z", num_timed, base_seed, camp)
        out[f"{p:g}"] = {"mean_ms": r["mean_ms"], "p999_ms": r["p999_ms"]}
        print(f"p={p:g}: mean={r['mean_ms']:.2f}ms p999={r['p999_ms']:.2f}ms")
    ref = {"5e-04": {"mean_ms": 3.55, "p999_ms": 272.5}, "0.001": {"mean_ms": 10.59, "p999_ms": 289.0}}
    verdict = {}
    for k, v in out.items():
        if k in ref:
            r_mean = v["mean_ms"] / ref[k]["mean_ms"]
            r_p999 = v["p999_ms"] / ref[k]["p999_ms"]
            verdict[k] = {
                "mean_ratio_vs_paper": r_mean,
                "p999_ratio_vs_paper": r_p999,
                "within_factor2": (
                    _within_factor_two_ratio(r_mean)
                    and _within_factor_two_ratio(r_p999)
                ),
            }
    camp.close({"gate": "A-timing", "measured": out, "published_ref": ref, "verdict": verdict})
    print("verdict:", json.dumps(verdict))


def cmd_gate_b(shots: int, base_seed: int = 20260829):
    camp = Campaign({
        "run": "gateB",
        "shots": shots,
        "base_seed": base_seed,
        "bposd_prior_model": BPOSD_PRIOR_MODEL,
    })
    _write_bposd_manifest(camp)
    r_bp = run_point("bp30+osd", 1e-3, "Z", shots, base_seed, camp)
    r_b8 = run_point("beam8", 1e-3, "Z", shots, base_seed, camp)
    ratio_point = r_b8["ler"] / r_bp["ler"] if r_bp["ler"] else float("inf")
    lo, hi = ratio_ci(r_b8["failures"], r_b8["shots"], r_bp["failures"], r_bp["shots"])
    if r_b8["failures"] == 0 and r_bp["failures"] == 0:
        # Zero failures in BOTH arms at this shot count: the ratio is
        # genuinely undefined and the data cannot distinguish the arms.
        # Per pre-statement semantics this is NOT a band verdict — record
        # an explicit under-determined outcome.
        verdict = {
            "beam8_vs_bposd_ler_ratio": None,
            "ratio_ci95": [lo, hi],
            "target_band": [0.87, 1.15],
            "in_band": None,
            "outcome": "UNDERDETERMINED_both_arms_zero_failures",
            "note": f"0 failures in {r_b8['shots']} shots on both arms; "
            "both decoders below resolvable LER at this budget. Increase "
            "shots or p; band check [0.87,1.15] requires resolvable arms.",
        }
    else:
        ratio_lo = r_b8["ler"] / max(r_bp["ler"], 1e-12)
        verdict = {
            "beam8_vs_bposd_ler_ratio": ratio_point,
            "ratio_ci95": [lo, hi],
            "target_band": [0.87, 1.15],
            "in_band": bool(0.87 <= ratio_lo <= 1.15 or (lo <= 1.15 and hi >= 0.87)),
        }
    camp.close({"gate": "B-beam8", "bposd": r_bp, "beam8": r_b8, "verdict": verdict})
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    import sys

    cmd = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if cmd == "smoke":
        cmd_smoke()
    elif cmd == "gateA_accuracy":
        shots = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000
        cmd_gate_a_accuracy(shots)
    elif cmd == "gateA_timing":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 10_000
        cmd_gate_a_timing(n)
    elif cmd == "gateB":
        shots = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000
        cmd_gate_b(shots)
    else:
        raise SystemExit(__doc__)
