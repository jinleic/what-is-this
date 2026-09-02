"""Resolve pinned Gate B at p=1e-3/Z with 1e6 shots per arm.

The two decoder arms preserve the frozen independent-sampling protocol from
``qldpc_dec.run_gate.cmd_gate_b``: each arm uses its existing
``sampling_seed(base_seed, p, basis, arm)`` stream. The beam8 arm is the
single-thread C++ port proven bit-exact to the pinned Python decoder on 1000
shared shots at both p=1e-3 and p=3e-3. BP+OSD decodes disjoint slices of one
frozen shot file in parallel; sharding changes execution only, not samples.

A campaign manifest and sample inventory are frozen before decoding. The run
closes only after every shard, artifact hash, shot count, and ratio verdict
validates.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
import re
import shutil
import struct
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import stim

TARGET = Path(__file__).resolve().parents[1]
SRC = TARGET / "src"
BEAM_DIR = SRC / "beam_cpp"
sys.path.insert(0, str(SRC))

from qldpc_dec.bootstrap import ratio_ci, wilson_interval  # noqa: E402
from qldpc_dec.bp_osd import BpOsdBatchDecoder  # noqa: E402
from qldpc_dec.circuits import circuit_sha, dem_sha, load_circuit, load_dem  # noqa: E402
from qldpc_dec.dem_matrices import dem_to_matrices  # noqa: E402
from qldpc_dec.runner import Campaign  # noqa: E402
from qldpc_dec.seeds import derive_seed, sampling_seed  # noqa: E402

P = 1e-3
BASIS = "Z"
ROUNDS = 12
BASE_SEED = 20260829
TARGET_BAND = (0.87, 1.15)
BOOTSTRAP_DRAWS = 100_000
BEAM_SOURCE = BEAM_DIR / "beam8.cpp"
BEAM_BINARY = BEAM_DIR / "beam8_cpp"
BOOTSTRAP_SOURCE = SRC / "qldpc_dec" / "bootstrap.py"
BP_OSD_SOURCE = SRC / "qldpc_dec" / "bp_osd.py"
DEM_MATRICES_SOURCE = SRC / "qldpc_dec" / "dem_matrices.py"
EQUIV_REPORTS = (
    BEAM_DIR / "scratch" / "equiv_p0.001_n1000_s20260830.json",
    BEAM_DIR / "scratch" / "equiv_p0.003_n1000_s20260830.json",
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while block := f.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def _u64_padded(packed: np.ndarray, nbits: int) -> np.ndarray:
    """Pad Stim's byte-packed rows to the u64-row format used by beam8_cpp."""
    packed = np.ascontiguousarray(packed, dtype=np.uint8)
    expected = (nbits + 7) // 8
    if packed.ndim != 2 or packed.shape[1] != expected:
        raise ValueError(
            f"packed width {packed.shape} does not match {nbits} bits ({expected} bytes)"
        )
    row_bytes = ((nbits + 63) // 64) * 8
    if row_bytes == expected:
        return packed
    out = np.zeros((len(packed), row_bytes), dtype=np.uint8)
    out[:, :expected] = packed
    return out


def _verify_bitpacked_sampling(circuit: stim.Circuit, seed: int) -> dict[str, object]:
    """Prove bit_packed=True preserves the seeded detector/observable stream."""
    n = 257
    det, obs = circuit.compile_detector_sampler(seed=seed).sample(
        n, separate_observables=True
    )
    det_p, obs_p = circuit.compile_detector_sampler(seed=seed).sample(
        n, separate_observables=True, bit_packed=True
    )
    det_u = np.unpackbits(det_p, axis=1, bitorder="little")[:, : det.shape[1]].astype(bool)
    obs_u = np.unpackbits(obs_p, axis=1, bitorder="little")[:, : obs.shape[1]].astype(bool)
    det_match = bool(np.array_equal(det, det_u))
    obs_match = bool(np.array_equal(obs, obs_u))
    if not det_match or not obs_match:
        raise AssertionError(
            f"Stim bit-packed sampler changed the seeded stream: det={det_match} obs={obs_match}"
        )
    return {"shots": n, "seed": seed, "detectors_match": True, "observables_match": True}


def _write_shots_file(
    path: Path, circuit: stim.Circuit, shots: int, seed: int
) -> dict[str, object]:
    sampler = circuit.compile_detector_sampler(seed=seed)
    det_p, obs_p = sampler.sample(
        shots, separate_observables=True, bit_packed=True
    )
    ndet = circuit.num_detectors
    nobs = circuit.num_observables
    det_rows = _u64_padded(det_p, ndet)
    obs_rows = _u64_padded(obs_p, nobs)
    with path.open("wb") as f:
        f.write(struct.pack("<QQQ", ndet, nobs, shots))
        f.write(det_rows.tobytes(order="C"))
        f.write(obs_rows.tobytes(order="C"))
    size = path.stat().st_size
    expected_size = 24 + shots * (det_rows.shape[1] + obs_rows.shape[1])
    if size != expected_size:
        raise AssertionError(f"shot file size {size} != expected {expected_size}")
    return {
        "path": path.name,
        "shots": shots,
        "sampling_seed": seed,
        "bytes": size,
        "sha256": sha256_file(path),
        "detector_row_bytes": int(det_rows.shape[1]),
        "observable_row_bytes": int(obs_rows.shape[1]),
    }


def _read_shots_header(path: Path) -> tuple[int, int, int, int, int]:
    with path.open("rb") as f:
        raw = f.read(24)
    if len(raw) != 24:
        raise ValueError(f"truncated shot file header: {path}")
    ndet, nobs, shots = struct.unpack("<QQQ", raw)
    det_bytes = ((ndet + 63) // 64) * 8
    obs_bytes = ((nobs + 63) // 64) * 8
    expected = 24 + shots * (det_bytes + obs_bytes)
    if path.stat().st_size != expected:
        raise ValueError(
            f"shot file size mismatch for {path}: {path.stat().st_size} != {expected}"
        )
    return ndet, nobs, shots, det_bytes, obs_bytes


def _load_shots_slice(path: Path, start: int, count: int) -> tuple[np.ndarray, np.ndarray]:
    ndet, nobs, shots, det_bytes, obs_bytes = _read_shots_header(path)
    if start < 0 or count < 0 or start + count > shots:
        raise ValueError(f"bad slice [{start}, {start + count}) for {shots} shots")
    det_p = np.memmap(
        path,
        dtype=np.uint8,
        mode="r",
        offset=24 + start * det_bytes,
        shape=(count, det_bytes),
    )
    obs_p = np.memmap(
        path,
        dtype=np.uint8,
        mode="r",
        offset=24 + shots * det_bytes + start * obs_bytes,
        shape=(count, obs_bytes),
    )
    det = np.unpackbits(det_p, axis=1, bitorder="little")[:, :ndet]
    obs = np.unpackbits(obs_p, axis=1, bitorder="little")[:, :nobs]
    return np.ascontiguousarray(det), np.ascontiguousarray(obs)


def _load_observables(path: Path) -> np.ndarray:
    _ndet, nobs, shots, det_bytes, obs_bytes = _read_shots_header(path)
    obs_p = np.memmap(
        path,
        dtype=np.uint8,
        mode="r",
        offset=24 + shots * det_bytes,
        shape=(shots, obs_bytes),
    )
    return np.unpackbits(obs_p, axis=1, bitorder="little")[:, :nobs]


def _write_failmask(path: Path, failed: np.ndarray) -> str:
    packed = np.packbits(np.asarray(failed, dtype=np.uint8), bitorder="little")
    path.write_bytes(packed.tobytes())
    return sha256_file(path)


def _decode_bposd_shard(
    campaign_dir_s: str, shard: int, start: int, count: int, shots_sha256: str
) -> dict[str, object]:
    campaign_dir = Path(campaign_dir_s)
    result_path = campaign_dir / f"bposd_shard_{shard:03d}.json"
    if result_path.exists():
        result = json.loads(result_path.read_text())
        expected = (shard, start, count, shots_sha256)
        got = (
            result.get("shard"),
            result.get("start"),
            result.get("shots"),
            result.get("source_shots_sha256"),
        )
        if got != expected:
            raise ValueError(f"stale shard result {result_path}: {got} != {expected}")
        return result

    shot_path = campaign_dir / "shots_bposd.bin"
    det, obs = _load_shots_slice(shot_path, start, count)
    dem = load_dem(P, BASIS)
    decoder = BpOsdBatchDecoder(dem)
    t0 = time.perf_counter()
    predictions = decoder.decode_batch(det)
    wall_s = time.perf_counter() - t0
    mismatch = np.asarray(predictions != obs)
    failed = mismatch.any(axis=1)
    mask_path = campaign_dir / f"bposd_failmask_{shard:03d}.bin"
    mask_sha = _write_failmask(mask_path, failed)
    result: dict[str, object] = {
        "arm": "bp30+osd",
        "shard": shard,
        "start": start,
        "shots": count,
        "failures": int(failed.sum()),
        "raw_logical_failure_rate": float(failed.mean()),
        "per_round_ler": float(failed.mean() / ROUNDS),
        "observable_mismatch_counts": mismatch.sum(axis=0).astype(int).tolist(),
        "wall_s": wall_s,
        "ms_per_shot": 1000.0 * wall_s / count,
        "source_shots_sha256": shots_sha256,
        "failure_mask": mask_path.name,
        "failure_mask_sha256": mask_sha,
    }
    result_path.write_text(json.dumps(result, indent=2))
    return result


def _parse_beam_predictions(path: Path, nobs: int, shots: int) -> np.ndarray:
    with path.open("rb") as f:
        header = f.read(12)
    if len(header) != 12:
        raise ValueError(f"truncated beam output header: {path}")
    magic, hdr_nobs, hdr_shots = struct.unpack("<III", header)
    if (magic, hdr_nobs, hdr_shots) != (0x31503842, nobs, shots):
        raise ValueError(
            f"beam header {(magic, hdr_nobs, hdr_shots)} != "
            f"{(0x31503842, nobs, shots)}"
        )
    row_bytes = ((nobs + 63) // 64) * 8
    expected = 12 + shots * row_bytes
    if path.stat().st_size != expected:
        raise ValueError(f"beam output size {path.stat().st_size} != {expected}")
    packed = np.memmap(path, dtype=np.uint8, mode="r", offset=12, shape=(shots, row_bytes))
    return np.unpackbits(packed, axis=1, bitorder="little")[:, :nobs]


def _beam_result(campaign_dir: Path, shots_sha256: str, wall_s: float) -> dict[str, object]:
    shot_path = campaign_dir / "shots_beam8.bin"
    pred_path = campaign_dir / "beam8_predictions.bin"
    _ndet, nobs, shots, _det_bytes, _obs_bytes = _read_shots_header(shot_path)
    predictions = _parse_beam_predictions(pred_path, nobs, shots)
    obs = _load_observables(shot_path)
    mismatch = np.asarray(predictions != obs)
    failed = mismatch.any(axis=1)
    mask_path = campaign_dir / "beam8_failmask.bin"
    mask_sha = _write_failmask(mask_path, failed)

    stderr_path = campaign_dir / "beam8.stderr.log"
    stderr = stderr_path.read_text() if stderr_path.exists() else ""
    match = re.search(r"=>\s*([0-9.eE+-]+)\s*ms/shot", stderr)
    decode_ms = float(match.group(1)) if match else 1000.0 * wall_s / shots
    return {
        "arm": "beam8-cpp-bit-exact",
        "shots": shots,
        "failures": int(failed.sum()),
        "raw_logical_failure_rate": float(failed.mean()),
        "per_round_ler": float(failed.mean() / ROUNDS),
        "observable_mismatch_counts": mismatch.sum(axis=0).astype(int).tolist(),
        "subprocess_wall_s": wall_s,
        "decode_ms_per_shot": decode_ms,
        "source_shots_sha256": shots_sha256,
        "predictions": pred_path.name,
        "predictions_sha256": sha256_file(pred_path),
        "failure_mask": mask_path.name,
        "failure_mask_sha256": mask_sha,
    }


def _aggregate_bposd(results: list[dict[str, object]], shots: int) -> dict[str, object]:
    ordered = sorted(results, key=lambda r: int(r["shard"]))
    if sum(int(r["shots"]) for r in ordered) != shots:
        raise AssertionError("BP+OSD shard shot counts do not sum to the campaign total")
    expected_start = 0
    for r in ordered:
        if int(r["start"]) != expected_start:
            raise AssertionError(f"BP+OSD shard gap/overlap at {expected_start}: {r}")
        expected_start += int(r["shots"])
    failures = sum(int(r["failures"]) for r in ordered)
    per_obs = np.sum(
        [np.asarray(r["observable_mismatch_counts"], dtype=np.int64) for r in ordered],
        axis=0,
    )
    raw = failures / shots
    core_s = sum(float(r["wall_s"]) for r in ordered)
    return {
        "arm": "bp30+osd",
        "shots": shots,
        "failures": failures,
        "raw_logical_failure_rate": raw,
        "raw_ler_ci95_wilson": list(wilson_interval(failures, shots)),
        "per_round_ler": raw / ROUNDS,
        "per_round_ler_ci95_wilson": [x / ROUNDS for x in wilson_interval(failures, shots)],
        "observable_mismatch_counts": per_obs.astype(int).tolist(),
        "decode_core_s": core_s,
        "decode_ms_per_shot": 1000.0 * core_s / shots,
        "shards": len(ordered),
    }


def _with_intervals(result: dict[str, object]) -> dict[str, object]:
    failures = int(result["failures"])
    shots = int(result["shots"])
    raw_ci = wilson_interval(failures, shots)
    result["raw_ler_ci95_wilson"] = list(raw_ci)
    result["per_round_ler_ci95_wilson"] = [x / ROUNDS for x in raw_ci]
    return result


def _ratio_verdict(
    beam: dict[str, object], bposd: dict[str, object], bootstrap_seed: int
) -> dict[str, object]:
    bf = int(beam["failures"])
    bn = int(beam["shots"])
    pf = int(bposd["failures"])
    pn = int(bposd["shots"])
    if bf == 0 or pf == 0:
        return {
            "outcome": "UNDERDETERMINED_zero_failures_in_at_least_one_arm",
            "beam8_over_bposd_ratio": None if pf == 0 else 0.0,
            "ratio_ci95": [0.0, None],
            "target_band": list(TARGET_BAND),
        }
    ratio = (bf / bn) / (pf / pn)
    lo, hi = ratio_ci(
        bf,
        bn,
        pf,
        pn,
        n_boot=BOOTSTRAP_DRAWS,
        rng=np.random.default_rng(bootstrap_seed),
    )
    point_in_band = bool(TARGET_BAND[0] <= ratio <= TARGET_BAND[1])
    ci_inside_band = bool(TARGET_BAND[0] <= lo and hi <= TARGET_BAND[1])
    ci_overlaps_band = bool(lo <= TARGET_BAND[1] and hi >= TARGET_BAND[0])
    if ci_inside_band:
        outcome = "REPRODUCED_ratio_CI_inside_band"
    elif not ci_overlaps_band:
        outcome = "NOT_REPRODUCED_ratio_CI_outside_band"
    else:
        outcome = "INCONCLUSIVE_ratio_CI_crosses_band"
    return {
        "outcome": outcome,
        "beam8_over_bposd_ratio": ratio,
        "ratio_ci95": [lo, float(hi) if np.isfinite(hi) else None],
        "ratio_ci95_upper_unbounded": not bool(np.isfinite(hi)),
        "target_band": list(TARGET_BAND),
        "point_in_band": point_in_band,
        "ci_inside_band": ci_inside_band,
        "ci_overlaps_band": ci_overlaps_band,
        "bootstrap": {
            "method": "independent-binomial percentile (frozen harness helper)",
            "draws": BOOTSTRAP_DRAWS,
            "seed": bootstrap_seed,
        },
    }


def _equivalence_evidence() -> list[dict[str, object]]:
    evidence: list[dict[str, object]] = []
    for path in EQUIV_REPORTS:
        report = json.loads(path.read_text())
        if report.get("mismatched_shots") != 0 or int(report.get("shots", 0)) < 1000:
            raise AssertionError(f"beam8 C++ equivalence prerequisite failed: {path}")
        evidence.append(
            {
                "path": str(path.relative_to(TARGET)),
                "sha256": sha256_file(path),
                "p": report["p"],
                "shots": report["shots"],
                "mismatched_shots": 0,
                "cpp_threads": report.get("cpp_threads", 1),
            }
        )
    return evidence


def _shard_bounds(shots: int, shards: int) -> list[tuple[int, int, int]]:
    if shots <= 0 or shards <= 0 or shards > shots:
        raise ValueError(f"invalid shots/shards: {shots}/{shards}")
    q, r = divmod(shots, shards)
    bounds: list[tuple[int, int, int]] = []
    start = 0
    for shard in range(shards):
        count = q + (1 if shard < r else 0)
        bounds.append((shard, start, count))
        start += count
    if start != shots:
        raise AssertionError((start, shots))
    return bounds


def run(args: argparse.Namespace) -> Path:
    if not BEAM_BINARY.exists() or not BEAM_SOURCE.exists():
        raise FileNotFoundError("beam8 source/binary missing; build and re-run equivalence first")
    evidence = _equivalence_evidence()
    shots = int(args.shots)
    shards = int(args.shards)
    workers = min(int(args.workers), shards)
    bounds = _shard_bounds(shots, shards)
    bposd_seed = sampling_seed(BASE_SEED, P, BASIS, "bp30+osd")
    beam_seed = sampling_seed(BASE_SEED, P, BASIS, "beam8")
    bootstrap_seed = derive_seed(BASE_SEED, "gateB-pinned-ratio-bootstrap", shots)
    config = {
        "run": "gateB-pinned-smoke" if args.smoke else "gateB-pinned-1e6",
        "p": P,
        "basis": BASIS,
        "rounds": ROUNDS,
        "shots_per_arm": shots,
        "base_seed": BASE_SEED,
        "sampler_seeds": {"bp30+osd": bposd_seed, "beam8": beam_seed},
        "sampling_design": "independent arms; exact existing sampling_seed streams",
        "bposd": {
            "max_iter": 30,
            "bp_method": "ms",
            "schedule": "parallel",
            "osd_method": "osd_cs",
            "osd_order": 10,
            "prior_model": "merged DEM heterogeneous error_channel",
            "shards": shards,
            "workers": workers,
        },
        "beam8": {
            "backend": "single-thread C++ bit-exact port",
            "beam_width": 8,
            "initial_iters": 30,
            "iters_per_round": 20,
            "max_rounds": 10,
            "num_results": 1,
        },
        "target_band_beam8_over_bposd": list(TARGET_BAND),
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "bootstrap_seed": bootstrap_seed,
        "ratio_bootstrap_zero_denominator": (
            "positive/zero retained as +inf; zero/zero undefined draws discarded"
        ),
    }
    root = TARGET / ("campaigns-smoke" if args.smoke else "campaigns")
    campaign = Campaign(config, root=root)
    campaign_dir = campaign.dir
    runner_copy = campaign_dir / "gateb_pinned.py"
    runner_copy.write_bytes(Path(__file__).read_bytes())
    bootstrap_copy = campaign_dir / "source_bootstrap.py"
    bootstrap_copy.write_bytes(BOOTSTRAP_SOURCE.read_bytes())
    bp_osd_copy = campaign_dir / "source_bp_osd.py"
    bp_osd_copy.write_bytes(BP_OSD_SOURCE.read_bytes())
    dem_matrices_copy = campaign_dir / "source_dem_matrices.py"
    dem_matrices_copy.write_bytes(DEM_MATRICES_SOURCE.read_bytes())

    circuit = load_circuit(P, BASIS)
    dem = load_dem(P, BASIS)
    if (dem.num_detectors, dem.num_errors, dem.num_observables) != (936, 8784, 12):
        raise AssertionError(
            f"wrong pinned DEM: {(dem.num_detectors, dem.num_errors, dem.num_observables)}"
        )
    (campaign_dir / "pinned.dem").write_text(str(dem))
    _H, _A, lam = dem_to_matrices(dem, merge=False)
    _merged_H, _merged_A, merged_priors = dem_to_matrices(dem, merge=True)
    prior_probe = BpOsdBatchDecoder(dem)
    configured_priors = np.asarray(prior_probe.decoder.channel_probs)
    np.testing.assert_array_equal(configured_priors, merged_priors)
    del prior_probe
    prior_check = {
        "columns": int(len(merged_priors)),
        "distinct_probabilities": int(np.unique(merged_priors).size),
        "minimum": float(merged_priors.min()),
        "maximum": float(merged_priors.max()),
        "sha256_float64": hashlib.sha256(
            np.asarray(merged_priors, dtype="<f8").tobytes()
        ).hexdigest(),
        "configured_channel_exact": True,
        "reference": {
            "package": "stimbposd==0.1.0",
            "commit": "7921f5eb1b358ff616f9822280c9961e83df06cb",
            "constructor_argument": "error_channel=list(priors)",
        },
    }
    np.save(campaign_dir / "pinned_lam.npy", lam)
    bitpack_check = _verify_bitpacked_sampling(circuit, bposd_seed)
    disk_free = shutil.disk_usage(campaign_dir).free
    if disk_free < 2_000_000_000:
        raise RuntimeError(f"need >=2 GB free for Gate B campaign; have {disk_free}")

    campaign.write_manifest(
        {
            "circuit_sha256": circuit_sha(P, BASIS),
            "dem_sha256": dem_sha(P, BASIS),
            "dem_dimensions": {
                "detectors": dem.num_detectors,
                "errors": dem.num_errors,
                "observables": dem.num_observables,
            },
            "runner_sha256": sha256_file(runner_copy),
            "bootstrap_source_sha256": sha256_file(bootstrap_copy),
            "bp_osd_source_sha256": sha256_file(bp_osd_copy),
            "dem_matrices_source_sha256": sha256_file(dem_matrices_copy),
            "bposd_prior_check": prior_check,
            "beam_source_sha256": sha256_file(BEAM_SOURCE),
            "beam_binary_sha256": sha256_file(BEAM_BINARY),
            "beam_equivalence_evidence": evidence,
            "versions": {
                "numpy": np.__version__,
                "stim": stim.__version__,
                "ldpc": importlib.metadata.version("ldpc"),
            },
            "bitpacked_sampler_check": bitpack_check,
            "pre_statement": "pre_statement.md Gate B revisions GB1-GB3",
        }
    )

    sample_inventory = {
        "sampling_frozen_before_decode": True,
        "arms": {
            "bp30+osd": _write_shots_file(
                campaign_dir / "shots_bposd.bin", circuit, shots, bposd_seed
            ),
            "beam8": _write_shots_file(
                campaign_dir / "shots_beam8.bin", circuit, shots, beam_seed
            ),
        },
    }
    (campaign_dir / "sample_inventory.json").write_text(
        json.dumps(sample_inventory, indent=2)
    )

    beam_stdout = (campaign_dir / "beam8.stdout.log").open("w")
    beam_stderr = (campaign_dir / "beam8.stderr.log").open("w")
    beam_cmd = [
        str(BEAM_BINARY),
        "--dem",
        str(campaign_dir / "pinned.dem"),
        "--lam",
        str(campaign_dir / "pinned_lam.npy"),
        "--shots",
        str(campaign_dir / "shots_beam8.bin"),
        "--out",
        str(campaign_dir / "beam8_predictions.bin"),
    ]
    beam_t0 = time.perf_counter()
    beam_proc = subprocess.Popen(beam_cmd, stdout=beam_stdout, stderr=beam_stderr)
    (campaign_dir / "launch.json").write_text(
        json.dumps(
            {
                "parent_pid": os.getpid(),
                "beam_pid": beam_proc.pid,
                "bposd_workers": workers,
                "bposd_shards": bounds,
                "beam_command": beam_cmd,
                "nice": os.nice(0),
            },
            indent=2,
        )
    )

    bposd_shots_sha = str(sample_inventory["arms"]["bp30+osd"]["sha256"])
    results: list[dict[str, object]] = []
    mp_context = multiprocessing.get_context("spawn")
    try:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers, mp_context=mp_context
        ) as pool:
            futures = {
                pool.submit(
                    _decode_bposd_shard,
                    str(campaign_dir),
                    shard,
                    start,
                    count,
                    bposd_shots_sha,
                ): shard
                for shard, start, count in bounds
            }
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                results.append(result)
                print(
                    f"BP+OSD shard {result['shard']}: "
                    f"{result['failures']}/{result['shots']} failures, "
                    f"{float(result['ms_per_shot']):.2f} ms/shot",
                    flush=True,
                )
        beam_returncode = beam_proc.wait()
    finally:
        beam_stdout.close()
        beam_stderr.close()
    beam_wall_s = time.perf_counter() - beam_t0
    if beam_returncode != 0:
        stderr_tail = (campaign_dir / "beam8.stderr.log").read_text()[-4000:]
        raise RuntimeError(f"beam8_cpp exited {beam_returncode}:\n{stderr_tail}")

    bposd = _aggregate_bposd(results, shots)
    beam = _with_intervals(
        _beam_result(
            campaign_dir,
            str(sample_inventory["arms"]["beam8"]["sha256"]),
            beam_wall_s,
        )
    )
    verdict = _ratio_verdict(beam, bposd, bootstrap_seed)
    campaign.append_result(bposd)
    campaign.append_result(beam)
    summary = {
        "gate": "B-beam8-pinned",
        "p": P,
        "basis": BASIS,
        "shots_per_arm": shots,
        "bposd": bposd,
        "beam8": beam,
        "verdict": verdict,
        "scope": (
            "beam8 first falsification target only; beam32/beam64 remain gated "
            "on this result per frozen pre_statement.md"
        ),
        "rate_semantics": {
            "raw_logical_failure_rate": "any-observable failure / shots",
            "per_round_ler": "raw_logical_failure_rate / 12",
            "ratio": "same under either convention; factor 1/12 cancels",
            "ratio_ci": (
                "independent-binomial percentile bootstrap; positive/zero "
                "draws retained as +inf"
            ),
        },
    }
    campaign.close(summary)
    print(json.dumps({"campaign": str(campaign_dir), **summary}, indent=2))
    return campaign_dir


def recover_unclosed(campaign_dir: Path) -> Path:
    """Validate and close a campaign after the known NumPy-bool JSON failure."""
    campaign_dir = campaign_dir.resolve()
    if (campaign_dir / "inventory.json").exists():
        raise ValueError(f"campaign is already closed: {campaign_dir}")
    launch = json.loads((campaign_dir / "launch.json").read_text())
    parent_pid = int(launch["parent_pid"])
    try:
        os.kill(parent_pid, 0)
    except ProcessLookupError:
        pass
    else:
        raise RuntimeError(f"campaign parent process {parent_pid} is still alive")

    manifest = json.loads((campaign_dir / "manifest.json").read_text())
    config = manifest["config"]
    shots = int(config["shots_per_arm"])
    with gzip.open(campaign_dir / "results.json.gz", "rt") as f:
        results = json.load(f)
    if len(results) != 2:
        raise ValueError(f"expected two completed arm results, got {len(results)}")
    bposd, beam = results
    if bposd["arm"] != "bp30+osd" or beam["arm"] != "beam8-cpp-bit-exact":
        raise ValueError("completed arm result order or labels are invalid")

    shard_count = int(config["bposd"]["shards"])
    shard_results = [
        json.loads((campaign_dir / f"bposd_shard_{i:03d}.json").read_text())
        for i in range(shard_count)
    ]
    if _aggregate_bposd(shard_results, shots) != bposd:
        raise ValueError("stored BP+OSD aggregate differs from validated shards")

    sample_inventory = json.loads((campaign_dir / "sample_inventory.json").read_text())
    for arm, filename in (
        ("bp30+osd", "shots_bposd.bin"),
        ("beam8", "shots_beam8.bin"),
    ):
        expected_sha = str(sample_inventory["arms"][arm]["sha256"])
        if sha256_file(campaign_dir / filename) != expected_sha:
            raise ValueError(f"{filename} differs from its frozen sample hash")

    predictions = _parse_beam_predictions(
        campaign_dir / str(beam["predictions"]), int(manifest["dem_dimensions"]["observables"]), shots
    )
    observables = _load_observables(campaign_dir / "shots_beam8.bin")
    mismatch = np.asarray(predictions != observables)
    failed = mismatch.any(axis=1)
    if int(failed.sum()) != int(beam["failures"]):
        raise ValueError("stored beam failure count differs from predictions")
    if mismatch.sum(axis=0).astype(int).tolist() != beam["observable_mismatch_counts"]:
        raise ValueError("stored beam observable counts differ from predictions")
    if sha256_file(campaign_dir / str(beam["predictions"])) != beam["predictions_sha256"]:
        raise ValueError("stored beam prediction hash differs from predictions file")
    if sha256_file(campaign_dir / str(beam["failure_mask"])) != beam["failure_mask_sha256"]:
        raise ValueError("stored beam failure-mask hash differs from failure-mask file")
    stored_mask = (campaign_dir / str(beam["failure_mask"])).read_bytes()
    if stored_mask != np.packbits(failed.astype(np.uint8), bitorder="little").tobytes():
        raise ValueError("stored beam failure mask differs from predictions")

    prior_model = config["bposd"].get("prior_model")
    if prior_model == "merged DEM heterogeneous error_channel":
        if not manifest["bposd_prior_check"]["configured_channel_exact"]:
            raise ValueError("corrected campaign lacks an exact configured-channel check")
        verdict = _ratio_verdict(beam, bposd, int(config["bootstrap_seed"]))
    else:
        verdict = {
            "outcome": "INVALIDATED_GB3_scalar_prior",
            "ratio_not_evaluated": True,
            "target_band": list(TARGET_BAND),
            "reason": (
                "BP+OSD used a scalar mean prior instead of the paper-pinned "
                "merged DEM heterogeneous error_channel"
            ),
        }

    recovery_runner = campaign_dir / "gateb_close_recovery.py"
    recovery_bootstrap = campaign_dir / "source_bootstrap_close_recovery.py"
    shutil.copy2(Path(__file__), recovery_runner)
    shutil.copy2(BOOTSTRAP_SOURCE, recovery_bootstrap)
    recovered_utc = datetime.now(timezone.utc).isoformat()
    recovery = {
        "recovered_utc": recovered_utc,
        "reason": (
            "initial close failed because NumPy bool values in the verdict "
            "were not JSON serializable"
        ),
        "runner_sha256": sha256_file(recovery_runner),
        "bootstrap_source_sha256": sha256_file(recovery_bootstrap),
        "validated_before_close": [
            "both frozen sample hashes",
            "all BP+OSD shard ranges and aggregate fields",
            "beam prediction, mismatch, and failure-mask hashes",
        ],
    }
    (campaign_dir / "recovery.json").write_text(json.dumps(recovery, indent=2))
    summary = {
        "gate": "B-beam8-pinned",
        "p": P,
        "basis": BASIS,
        "shots_per_arm": shots,
        "bposd": bposd,
        "beam8": beam,
        "verdict": verdict,
        "scope": (
            "beam8 first falsification target only; beam32/beam64 remain gated "
            "on this result per frozen pre_statement.md"
        ),
        "rate_semantics": {
            "raw_logical_failure_rate": "any-observable failure / shots",
            "per_round_ler": "raw_logical_failure_rate / 12",
            "ratio": "same under either convention; factor 1/12 cancels",
            "ratio_ci": (
                "independent-binomial percentile bootstrap; positive/zero "
                "draws retained as +inf"
            ),
        },
        "closure_recovery": recovery,
    }
    (campaign_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    inventory = {
        "files": sorted(
            p.name for p in campaign_dir.iterdir() if p.name != "inventory.json"
        )
        + ["inventory.json"],
        "closed_utc": recovered_utc,
        "results_count": len(results),
        "policy": "immutable-after-close",
    }
    (campaign_dir / "inventory.json").write_text(json.dumps(inventory, indent=2))
    print(json.dumps({"campaign": str(campaign_dir), **summary}, indent=2))
    return campaign_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shots", type=int, default=1_000_000)
    parser.add_argument("--shards", type=int, default=20)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--recover", type=Path)
    args = parser.parse_args()
    if args.recover is not None:
        recover_unclosed(args.recover)
        return
    if args.smoke and args.shots == 1_000_000:
        args.shots = 200
    if args.smoke and args.shards == 20:
        args.shards = 2
    if args.smoke and args.workers == 16:
        args.workers = 2
    try:
        os.nice(10)
    except OSError:
        pass
    run(args)


if __name__ == "__main__":
    main()
