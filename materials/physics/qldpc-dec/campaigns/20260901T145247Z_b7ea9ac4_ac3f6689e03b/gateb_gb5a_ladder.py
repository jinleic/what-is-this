"""GB5/GB5a Gate-B ladder + paired width instrument (pre_statement.md GB1-GB5a).

Two instruments, one campaign, at p=1e-3 basis Z:

1. LADDER (frozen bands, independent arms). Each rung decodes its own frozen
   stream sampling_seed(20260829, 1e-3, "Z", rung); the BP+OSD denominator is
   decoded once on sampling_seed(..., "bp30+osd") and shared by all three
   rungs (disclosed dependence). Ratio := beam raw failure rate / BP+OSD raw
   failure rate; bands in that convention:
       beam8  [0.87, 1.15]        (published "same LER")
       beam32 [1/7, 1/4.5]        (published "5.6x lower")
       beam64 [1/22, 1/11]        (published "17x lower")
   Three-way rule (GB1) on the GB2 sparse-tail bootstrap: CI wholly inside ->
   REPRODUCED; wholly outside -> NOT_REPRODUCED; overlapping -> INCONCLUSIVE.

2. PAIRED INSTRUMENT (GB5a). All three beam parameter sets additionally decode
   the SAME shared stream that BP+OSD decoded. On identical shots this gives
   per-decoder failure sets, pairwise discordance counts n(A fails, B ok),
   exact McNemar p-values, and prediction-identity hashes. The published
   ladder requires beam32 to fix ~82% and beam64 ~94% of beam8's failures;
   zero discordance refutes that as a sample identity, not a CI statement.

Launch gate: gb5eq_*.json equivalence reports (Python vs the parameter-wired
C++ binary) must exist with zero mismatches for every rung at both p=1e-3 and
p=3e-3, and their hashes are frozen into the manifest.

Every decode step is idempotent and hash-checked within a run. A campaign
directory is created fresh per invocation (Campaign contract), so an aborted
run restarts from new samples; the frozen seeds make that deterministic, only
costly.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import json
import math
import multiprocessing
import os
import re
import struct
import subprocess
import sys
import time
from pathlib import Path
from types import MappingProxyType

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
BOOTSTRAP_DRAWS = 100_000
CHUNK = 1_000_000

BEAM_SOURCE = BEAM_DIR / "beam8.cpp"
BEAM_BINARY = BEAM_DIR / "beam8_cpp"
BOOTSTRAP_SOURCE = SRC / "qldpc_dec" / "bootstrap.py"
BP_OSD_SOURCE = SRC / "qldpc_dec" / "bp_osd.py"
DEM_MATRICES_SOURCE = SRC / "qldpc_dec" / "dem_matrices.py"
BEAM_SEARCH_SOURCE = SRC / "qldpc_dec" / "beam_search.py"

RUNGS = MappingProxyType({
    "beam8": MappingProxyType({
        "band": (0.87, 1.15),
        "published": "same LER as bp30+osd",
        "beam_width": 8, "initial_iters": 30, "iters_per_round": 20, "max_rounds": 10,
    }),
    "beam32": MappingProxyType({
        "band": (1.0 / 7.0, 1.0 / 4.5),
        "published": "5.6x lower LER than bp30+osd",
        "beam_width": 32, "initial_iters": 40, "iters_per_round": 30, "max_rounds": 10,
    }),
    "beam64": MappingProxyType({
        "band": (1.0 / 22.0, 1.0 / 11.0),
        "published": "17x lower LER than bp30+osd",
        "beam_width": 64, "initial_iters": 40, "iters_per_round": 30, "max_rounds": 20,
    }),
})
IONQ_BPOSD = MappingProxyType({
    "max_iter": 30, "bp_method": "ms", "schedule": "parallel",
    "osd_method": "osd_cs", "osd_order": 10,
})
EQUIV_MATRIX = tuple((rung, p, n) for rung in RUNGS for p, n in ((1e-3, 2000), (3e-3, 300)))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _u64_padded(packed: np.ndarray, nbits: int) -> np.ndarray:
    packed = np.ascontiguousarray(packed, dtype=np.uint8)
    expected = (nbits + 7) // 8
    if packed.ndim != 2 or packed.shape[1] != expected:
        raise ValueError(f"packed width {packed.shape} does not match {nbits} bits")
    row_bytes = ((nbits + 63) // 64) * 8
    if row_bytes == expected:
        return packed
    out = np.zeros((len(packed), row_bytes), dtype=np.uint8)
    out[:, :expected] = packed
    return out


def write_shots(path: Path, circuit: stim.Circuit, shots: int, seed: int) -> dict[str, object]:
    """Sample once from a freshly seeded sampler with bounded memory."""
    ndet, nobs = circuit.num_detectors, circuit.num_observables
    det_row = ((ndet + 63) // 64) * 8
    obs_row = ((nobs + 63) // 64) * 8
    expected = 24 + shots * (det_row + obs_row)
    sampler = circuit.compile_detector_sampler(seed=seed)
    with path.open("w+b") as f:
        f.write(struct.pack("<QQQ", ndet, nobs, shots))
        f.truncate(expected)
        start = 0
        while start < shots:
            take = min(CHUNK, shots - start)
            det_p, obs_p = sampler.sample(
                take, separate_observables=True, bit_packed=True)
            det = _u64_padded(det_p, ndet)
            obs = _u64_padded(obs_p, nobs)
            f.seek(24 + start * det_row)
            f.write(det.tobytes(order="C"))
            f.seek(24 + shots * det_row + start * obs_row)
            f.write(obs.tobytes(order="C"))
            start += take
        f.flush()
        os.fsync(f.fileno())
    size = path.stat().st_size
    if size != expected:
        raise AssertionError(f"shot file size {size} != expected {expected}")
    return {
        "path": path.name, "shots": shots, "sampling_seed": seed,
        "bytes": size, "sha256": sha256_file(path),
        "detector_row_bytes": det_row, "observable_row_bytes": obs_row,
    }


def read_shots_header(path: Path) -> tuple[int, int, int, int, int]:
    with path.open("rb") as f:
        raw = f.read(24)
    if len(raw) != 24:
        raise ValueError(f"truncated shot file header: {path}")
    ndet, nobs, shots = struct.unpack("<QQQ", raw)
    det_bytes = ((ndet + 63) // 64) * 8
    obs_bytes = ((nobs + 63) // 64) * 8
    if path.stat().st_size != 24 + shots * (det_bytes + obs_bytes):
        raise ValueError(f"shot file size mismatch for {path}")
    return ndet, nobs, shots, det_bytes, obs_bytes


def load_shots_slice(path: Path, start: int, count: int) -> tuple[np.ndarray, np.ndarray]:
    ndet, nobs, shots, det_bytes, obs_bytes = read_shots_header(path)
    if start < 0 or count < 0 or start + count > shots:
        raise ValueError(f"bad slice [{start}, {start + count}) for {shots} shots")
    det_p = np.memmap(path, dtype=np.uint8, mode="r",
                      offset=24 + start * det_bytes, shape=(count, det_bytes))
    obs_p = np.memmap(path, dtype=np.uint8, mode="r",
                      offset=24 + shots * det_bytes + start * obs_bytes,
                      shape=(count, obs_bytes))
    det = np.unpackbits(det_p, axis=1, bitorder="little")[:, :ndet]
    obs = np.unpackbits(obs_p, axis=1, bitorder="little")[:, :nobs]
    return np.ascontiguousarray(det), np.ascontiguousarray(obs)


def load_observables_chunk(path: Path, start: int, count: int) -> np.ndarray:
    ndet, nobs, shots, det_bytes, obs_bytes = read_shots_header(path)
    obs_p = np.memmap(path, dtype=np.uint8, mode="r",
                      offset=24 + shots * det_bytes + start * obs_bytes,
                      shape=(count, obs_bytes))
    return np.unpackbits(obs_p, axis=1, bitorder="little")[:, :nobs]


def load_predictions_chunk(path: Path, nobs: int, shots: int, start: int, count: int) -> np.ndarray:
    row_bytes = ((nobs + 63) // 64) * 8
    packed = np.memmap(path, dtype=np.uint8, mode="r",
                       offset=12 + start * row_bytes, shape=(count, row_bytes))
    return np.unpackbits(packed, axis=1, bitorder="little")[:, :nobs]


def check_predictions_header(path: Path, nobs: int, shots: int) -> None:
    with path.open("rb") as f:
        header = f.read(12)
    magic, hdr_nobs, hdr_shots = struct.unpack("<III", header)
    if (magic, hdr_nobs, hdr_shots) != (0x31503842, nobs, shots):
        raise ValueError(f"beam header {(magic, hdr_nobs, hdr_shots)} mismatch for {path}")
    row_bytes = ((nobs + 63) // 64) * 8
    if path.stat().st_size != 12 + shots * row_bytes:
        raise ValueError(f"beam output size mismatch for {path}")


def failure_mask_from_predictions(
    shot_path: Path, pred_path: Path
) -> tuple[np.ndarray, np.ndarray]:
    """Chunked per-shot failure mask + per-observable mismatch counts."""
    _ndet, nobs, shots, _db, _ob = read_shots_header(shot_path)
    check_predictions_header(pred_path, nobs, shots)
    failed = np.zeros(shots, dtype=bool)
    per_obs = np.zeros(nobs, dtype=np.int64)
    for start in range(0, shots, CHUNK):
        count = min(CHUNK, shots - start)
        preds = load_predictions_chunk(pred_path, nobs, shots, start, count)
        obs = load_observables_chunk(shot_path, start, count)
        mism = preds != obs
        failed[start:start + count] = mism.any(axis=1)
        per_obs += mism.sum(axis=0).astype(np.int64)
    return failed, per_obs


def write_failmask(path: Path, failed: np.ndarray) -> str:
    path.write_bytes(np.packbits(failed.astype(np.uint8), bitorder="little").tobytes())
    return sha256_file(path)


def read_failmask(path: Path, shots: int) -> np.ndarray:
    raw = np.frombuffer(path.read_bytes(), dtype=np.uint8)
    return np.unpackbits(raw, bitorder="little")[:shots].astype(bool)


def shard_bounds(shots: int, shards: int) -> list[tuple[int, int, int]]:
    if shots <= 0 or shards <= 0 or shards > shots:
        raise ValueError(f"invalid shots/shards: {shots}/{shards}")
    q, r = divmod(shots, shards)
    out, start = [], 0
    for shard in range(shards):
        count = q + (1 if shard < r else 0)
        out.append((shard, start, count))
        start += count
    if start != shots:
        raise AssertionError((start, shots))
    return out


def decode_bposd_shard(
    campaign_dir_s: str, shard: int, start: int, count: int, shots_sha256: str
) -> dict[str, object]:
    campaign_dir = Path(campaign_dir_s)
    result_path = campaign_dir / f"bposd_shard_{shard:04d}.json"
    if result_path.exists():
        result = json.loads(result_path.read_text())
        got = (result.get("shard"), result.get("start"), result.get("shots"),
               result.get("source_shots_sha256"))
        if got != (shard, start, count, shots_sha256):
            raise ValueError(f"stale shard result {result_path}: {got}")
        return result
    det, obs = load_shots_slice(campaign_dir / "shots_bposd.bin", start, count)
    dem = load_dem(P, BASIS)
    decoder = BpOsdBatchDecoder(dem)
    t0 = time.perf_counter()
    predictions = decoder.decode_batch(det)
    wall_s = time.perf_counter() - t0
    mismatch = np.asarray(predictions != obs)
    failed = mismatch.any(axis=1)
    mask_path = campaign_dir / f"bposd_failmask_{shard:04d}.bin"
    mask_sha = write_failmask(mask_path, failed)
    result: dict[str, object] = {
        "arm": "bp30+osd", "shard": shard, "start": start, "shots": count,
        "failures": int(failed.sum()),
        "observable_mismatch_counts": mismatch.sum(axis=0).astype(int).tolist(),
        "wall_s": wall_s, "ms_per_shot": 1000.0 * wall_s / count,
        "source_shots_sha256": shots_sha256,
        "failure_mask": mask_path.name, "failure_mask_sha256": mask_sha,
    }
    result_path.write_text(json.dumps(result, indent=2))
    return result


def run_beam(
    campaign_dir: Path, rung: str, shots_file: str, out_name: str, log_tag: str,
    threads: int, source_shots_sha256: str,
) -> dict[str, object]:
    """Decode one shot file with one rung's pinned parameters (idempotent)."""
    params = RUNGS[rung]
    expected_cfg = (
        f"config: beam_width={params['beam_width']} "
        f"initial_iters={params['initial_iters']} "
        f"iters_per_round={params['iters_per_round']} "
        f"max_rounds={params['max_rounds']} num_results=1"
    )
    binary_sha = sha256_file(BEAM_BINARY)
    pred_path = campaign_dir / out_name
    meta_path = campaign_dir / f"{out_name}.meta.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        valid = (
            pred_path.is_file()
            and meta.get("predictions_sha256") == sha256_file(pred_path)
            and meta.get("beam_binary_sha256") == binary_sha
            and meta.get("rung") == rung
            and meta.get("source_shots_sha256") == source_shots_sha256
            and meta.get("shots_file") == shots_file
            and meta.get("config_line") == expected_cfg
            and meta.get("threads") == threads
        )
        if valid:
            return meta
        raise ValueError(f"stale beam predictions or metadata for {out_name}")
    cmd = [
        str(BEAM_BINARY),
        "--dem", str(campaign_dir / "pinned.dem"),
        "--lam", str(campaign_dir / "pinned_lam.npy"),
        "--shots", str(campaign_dir / shots_file),
        "--out", str(pred_path),
        f"--beam-width={params['beam_width']}",
        f"--initial-iters={params['initial_iters']}",
        f"--iters-per-round={params['iters_per_round']}",
        f"--max-rounds={params['max_rounds']}",
        "--threads", str(threads),
    ]
    stdout_path = campaign_dir / f"{log_tag}.stdout.log"
    stderr_path = campaign_dir / f"{log_tag}.stderr.log"
    t0 = time.perf_counter()
    with stdout_path.open("w") as so, stderr_path.open("w") as se:
        rc = subprocess.call(cmd, stdout=so, stderr=se)
    wall = time.perf_counter() - t0
    if rc != 0:
        raise RuntimeError(f"{rung} decode failed ({rc}):\n{stderr_path.read_text()[-4000:]}")
    stderr = stderr_path.read_text()
    cfg = next((l.strip() for l in stderr.splitlines() if l.startswith("config:")), "")
    if cfg != expected_cfg:
        raise AssertionError(f"binary ran the wrong configuration:\n  got {cfg}\n  want {expected_cfg}")
    stats = next((l.strip() for l in stderr.splitlines() if l.startswith("stats:")), "")
    ms_match = re.search(r"=>\s*([0-9.eE+-]+)\s*ms/shot", stderr)
    meta = {
        "rung": rung,
        "shots_file": shots_file,
        "predictions": out_name,
        "predictions_sha256": sha256_file(pred_path),
        "beam_binary_sha256": binary_sha,
        "source_shots_sha256": source_shots_sha256,
        "config_line": cfg,
        "stats_line": stats,
        "threads": threads,
        "wall_s": wall,
        "ms_per_shot_reported": float(ms_match.group(1)) if ms_match else None,
        "command": cmd,
    }
    meta_path.write_text(json.dumps(meta, indent=2))
    return meta


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p-value for discordant pairs (b, c)."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2 ** n)
    return float(min(1.0, 2.0 * tail))


def equivalence_evidence() -> list[dict[str, object]]:
    """Launch gate: zero mismatches for every rung at both noise points."""
    out: list[dict[str, object]] = []
    for rung, p, n in EQUIV_MATRIX:
        path = BEAM_DIR / "evidence" / f"gb5eq_{rung}_p{p:g}_n{n}.json"
        if not path.exists():
            raise FileNotFoundError(f"missing equivalence report {path}")
        report = json.loads(path.read_text())
        identity = (
            report.get("rung") == rung
            and float(report.get("p", -1.0)) == p
            and report.get("basis") == BASIS
            and int(report.get("shots", -1)) == n
            and report.get("circuit_sha256") == circuit_sha(p, BASIS)
            and report.get("dem_sha256") == dem_sha(p, BASIS)
            and report.get("beam_source_sha256") == sha256_file(BEAM_SOURCE)
        )
        if not identity:
            raise AssertionError(f"equivalence report identity mismatch for {path}")
        if report["mismatched_shots"] != 0:
            raise AssertionError(f"equivalence FAILED for {path}: {report['mismatched_shots']}")
        if report["beam_binary_sha256"] != sha256_file(BEAM_BINARY):
            raise AssertionError(f"equivalence report {path} is for a different binary")
        pinned = {
            **{k: RUNGS[rung][k] for k in
               ("beam_width", "initial_iters", "iters_per_round", "max_rounds")},
            "num_results": 1,
        }
        if {k: report["parameters"].get(k) for k in pinned} != pinned:
            raise AssertionError(
                f"equivalence report {path} validated {report['parameters']}, "
                f"not the pinned {rung} parameters {pinned}")
        out.append({
            "path": str(path.relative_to(TARGET)), "sha256": sha256_file(path),
            "rung": rung, "p": p, "shots": report["shots"],
            "mismatched_shots": 0, "parameters": report["parameters"],
            "cpp_config_line": report["cpp_config_line"],
        })
    return out


def arm_summary(rung: str, shots: int, failures: int, per_obs: list[int],
                meta: dict[str, object], mask_name: str, mask_sha: str) -> dict[str, object]:
    raw = failures / shots
    ci = wilson_interval(failures, shots)
    return {
        "arm": f"{rung}-cpp-bit-exact", "rung": rung, "shots": shots, "failures": failures,
        "raw_logical_failure_rate": raw, "raw_ler_ci95_wilson": list(ci),
        "per_round_ler": raw / ROUNDS,
        "per_round_ler_ci95_wilson": [x / ROUNDS for x in ci],
        "observable_mismatch_counts": per_obs,
        "decoder_params": {k: RUNGS[rung][k] for k in
                           ("beam_width", "initial_iters", "iters_per_round", "max_rounds")},
        "decode": meta,
        "failure_mask": mask_name, "failure_mask_sha256": mask_sha,
    }


def ratio_verdict(beam: dict[str, object], bposd: dict[str, object],
                  band: tuple[float, float], bootstrap_seed: int) -> dict[str, object]:
    bf, bn = int(beam["failures"]), int(beam["shots"])
    pf, pn = int(bposd["failures"]), int(bposd["shots"])
    base = {"target_band": [float(band[0]), float(band[1])],
            "bootstrap": {"method": "independent-binomial percentile (frozen harness helper)",
                          "draws": BOOTSTRAP_DRAWS, "seed": bootstrap_seed}}
    if bf == 0 or pf == 0:
        # Parametric-bootstrap draws from p_hat = 0 are degenerate (every draw
        # is 0), so a CI computed here would manufacture a false decisive
        # verdict. Keep the frozen GB1 UNDERDETERMINED outcome and report a
        # non-degenerate Wilson-based one-sided bound alongside it.
        num_hi = wilson_interval(bf, bn)[1]
        den_lo = wilson_interval(pf, pn)[0]
        return {**base, "outcome": "UNDERDETERMINED_zero_failures_in_at_least_one_arm",
                "beam_over_bposd_ratio": None, "ratio_ci95": [0.0, None],
                "point_in_band": False, "ci_inside_band": False, "ci_overlaps_band": True,
                "wilson_ratio_upper_bound": (float(num_hi / den_lo) if den_lo > 0 else None),
                "note": "zero failures in one arm; parametric bootstrap is degenerate here"}
    ratio = (bf / bn) / (pf / pn)
    lo, hi = ratio_ci(bf, bn, pf, pn, n_boot=BOOTSTRAP_DRAWS,
                      rng=np.random.default_rng(bootstrap_seed))
    inside = bool(band[0] <= lo and hi <= band[1])
    overlaps = bool(lo <= band[1] and hi >= band[0])
    outcome = ("REPRODUCED_ratio_CI_inside_band" if inside else
               "INCONCLUSIVE_ratio_CI_crosses_band" if overlaps else
               "NOT_REPRODUCED_ratio_CI_outside_band")
    return {**base, "outcome": outcome, "beam_over_bposd_ratio": float(ratio),
            "ratio_ci95": [float(lo), float(hi) if np.isfinite(hi) else None],
            "ratio_ci95_upper_unbounded": not bool(np.isfinite(hi)),
            "point_in_band": bool(band[0] <= ratio <= band[1]),
            "ci_inside_band": inside, "ci_overlaps_band": overlaps}


def ladder_verdict(rung_verdicts: dict[str, dict[str, object]]) -> dict[str, object]:
    outcomes = {r: str(v["outcome"]) for r, v in rung_verdicts.items()}
    if any(o.startswith("NOT_REPRODUCED") for o in outcomes.values()):
        ladder = "REFUTED_at_least_one_published_band_excluded"
    elif all(o.startswith("REPRODUCED") for o in outcomes.values()):
        ladder = "REPRODUCED_all_published_bands_recovered"
    else:
        ladder = "INCONCLUSIVE_no_band_excluded_or_recovered_for_every_rung"
    return {"ladder_outcome": ladder, "rung_outcomes": outcomes}


def run(args: argparse.Namespace) -> Path:
    shots = int(args.shots)
    shards = int(args.shards)
    workers = int(args.workers)
    threads = int(args.threads)
    evidence = equivalence_evidence()

    config = {
        "run": "gateB-gb5a-ladder" + ("-smoke" if args.smoke else ""),
        "p": P, "basis": BASIS, "rounds": ROUNDS, "base_seed": BASE_SEED,
        "shots_per_arm": shots,
        "rungs": {r: {"band": [float(RUNGS[r]["band"][0]), float(RUNGS[r]["band"][1])],
                      "published": RUNGS[r]["published"],
                      **{k: RUNGS[r][k] for k in ("beam_width", "initial_iters",
                                                  "iters_per_round", "max_rounds")},
                      "num_results": 1}
                  for r in RUNGS},
        "bposd": {**dict(IONQ_BPOSD),
                  "prior_model": "merged DEM heterogeneous error_channel",
                  "shards": shards, "workers": workers,
                  "denominator_shared_by_all_rungs": True},
        "paired_instrument": {
            "shared_stream": "bp30+osd",
            "decoders": ["bp30+osd", *RUNGS],
            "reports": ["failure sets", "pairwise discordance", "exact McNemar", "identity hashes"],
        },
        "bootstrap": {"draws": BOOTSTRAP_DRAWS,
                      "zero_denominator": "positive/zero retained as +inf; zero/zero discarded"},
        "beam_threads": threads,
        "pre_statement": "GB1-GB5a",
    }
    root = TARGET / ("campaigns-smoke" if args.smoke else "campaigns")
    campaign = Campaign(config, root=root)
    cdir = campaign.dir
    print(f"campaign: {cdir}", flush=True)
    for src, dst in ((Path(__file__), "gateb_ladder.py"),
                     (BOOTSTRAP_SOURCE, "source_bootstrap.py"),
                     (BP_OSD_SOURCE, "source_bp_osd.py"),
                     (DEM_MATRICES_SOURCE, "source_dem_matrices.py"),
                     (BEAM_SEARCH_SOURCE, "source_beam_search.py"),
                     (BEAM_SOURCE, "source_beam8.cpp"),
                     (TARGET / "pre_statement.md", "pre_statement.md")):
        (cdir / dst).write_bytes(src.read_bytes())

    circuit = load_circuit(P, BASIS)
    dem = load_dem(P, BASIS)
    if (dem.num_detectors, dem.num_errors, dem.num_observables) != (936, 8784, 12):
        raise AssertionError("wrong pinned DEM shape")
    (cdir / "pinned.dem").write_text(str(dem))
    _H, _A, lam = dem_to_matrices(dem, merge=False)
    np.save(cdir / "pinned_lam.npy", lam)
    _mH, _mA, merged_priors = dem_to_matrices(dem, merge=True)
    probe = BpOsdBatchDecoder(dem)
    np.testing.assert_array_equal(np.asarray(probe.decoder.channel_probs), merged_priors)
    del probe
    prior_check = {
        "columns": int(len(merged_priors)),
        "distinct_probabilities": int(np.unique(merged_priors).size),
        "minimum": float(merged_priors.min()), "maximum": float(merged_priors.max()),
        "sha256_float64": hashlib.sha256(np.asarray(merged_priors, dtype="<f8").tobytes()).hexdigest(),
        "configured_channel_exact": True,
        "reference": {"package": "stimbposd==0.1.0",
                      "commit": "7921f5eb1b358ff616f9822280c9961e83df06cb",
                      "constructor_argument": "error_channel=list(priors)"},
    }

    print("sampling frozen streams ...", flush=True)
    inventory: dict[str, object] = {"sampling_frozen_before_decode": True, "arms": {}}
    arms: dict[str, object] = inventory["arms"]  # type: ignore[assignment]
    bp_seed = sampling_seed(BASE_SEED, P, BASIS, "bp30+osd")
    arms["bp30+osd"] = write_shots(cdir / "shots_bposd.bin", circuit, shots, bp_seed)
    for rung in RUNGS:
        seed = sampling_seed(BASE_SEED, P, BASIS, rung)
        arms[rung] = write_shots(cdir / f"shots_{rung}.bin", circuit, shots, seed)
    (cdir / "sample_inventory.json").write_text(json.dumps(inventory, indent=2))

    campaign.write_manifest({
        "circuit_sha256": circuit_sha(P, BASIS), "dem_sha256": dem_sha(P, BASIS),
        "dem_dimensions": {"detectors": dem.num_detectors, "errors": dem.num_errors,
                           "observables": dem.num_observables},
        "runner_sha256": sha256_file(cdir / "gateb_ladder.py"),
        "bootstrap_source_sha256": sha256_file(cdir / "source_bootstrap.py"),
        "bp_osd_source_sha256": sha256_file(cdir / "source_bp_osd.py"),
        "dem_matrices_source_sha256": sha256_file(cdir / "source_dem_matrices.py"),
        "beam_search_source_sha256": sha256_file(cdir / "source_beam_search.py"),
        "beam_source_sha256": sha256_file(BEAM_SOURCE),
        "beam_binary_sha256": sha256_file(BEAM_BINARY),
        "beam_equivalence_evidence": evidence,
        "bposd_prior_check": prior_check,
        "versions": {"numpy": np.__version__, "stim": stim.__version__,
                     "ldpc": importlib.metadata.version("ldpc")},
        "pre_statement": "pre_statement.md Gate B revisions GB1-GB5a",
    })

    # ---- BP+OSD denominator (shared by all rungs) ----
    print(f"BP+OSD: {shots} shots over {shards} shards on {workers} workers ...", flush=True)
    bp_sha = str(arms["bp30+osd"]["sha256"])  # type: ignore[index]
    bounds = shard_bounds(shots, shards)
    results: list[dict[str, object]] = []
    ctx = multiprocessing.get_context("spawn")
    t_bp = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as pool:
        futs = {pool.submit(decode_bposd_shard, str(cdir), s, st, c, bp_sha): s
                for s, st, c in bounds}
        done = 0
        for fut in concurrent.futures.as_completed(futs):
            res = fut.result()
            results.append(res)
            done += 1
            if done % 10 == 0 or done == len(bounds):
                tot = sum(int(r["failures"]) for r in results)
                print(f"  shards {done}/{len(bounds)} | failures so far {tot}", flush=True)
    bp_wall = time.perf_counter() - t_bp
    results.sort(key=lambda r: int(r["shard"]))
    if sum(int(r["shots"]) for r in results) != shots:
        raise AssertionError("BP+OSD shards do not sum to the campaign shot count")
    expected = 0
    for r in results:
        if int(r["start"]) != expected:
            raise AssertionError(f"BP+OSD shard gap at {expected}")
        expected += int(r["shots"])
    bp_failures = sum(int(r["failures"]) for r in results)
    for r in results:
        mask_path = cdir / str(r["failure_mask"])
        if sha256_file(mask_path) != str(r["failure_mask_sha256"]):
            raise AssertionError(f"shard failure mask changed on disk: {mask_path}")
    bp_mask = np.concatenate([
        read_failmask(cdir / str(r["failure_mask"]), int(r["shots"])) for r in results])
    if bp_mask.size != shots or int(bp_mask.sum()) != bp_failures:
        raise AssertionError("stitched BP+OSD mask disagrees with shard counts")
    bp_mask_sha = write_failmask(cdir / "bposd_failmask.bin", bp_mask)
    per_obs = np.sum([np.asarray(r["observable_mismatch_counts"], dtype=np.int64)
                      for r in results], axis=0)
    raw = bp_failures / shots
    ci = wilson_interval(bp_failures, shots)
    core_s = sum(float(r["wall_s"]) for r in results)
    bposd = {
        "arm": "bp30+osd", "shots": shots, "failures": bp_failures,
        "raw_logical_failure_rate": raw, "raw_ler_ci95_wilson": list(ci),
        "per_round_ler": raw / ROUNDS,
        "per_round_ler_ci95_wilson": [x / ROUNDS for x in ci],
        "observable_mismatch_counts": per_obs.astype(int).tolist(),
        "decode_core_s": core_s, "decode_ms_per_shot": 1000.0 * core_s / shots,
        "wall_s": bp_wall, "shards": len(results), "workers": workers,
        "decoder": dict(IONQ_BPOSD),
        "failure_mask": "bposd_failmask.bin", "failure_mask_sha256": bp_mask_sha,
        "source_shots_sha256": bp_sha,
    }
    campaign.append_result(bposd)
    print(f"BP+OSD: {bp_failures}/{shots} failures "
          f"({1000.0 * core_s / shots:.2f} ms/shot core, {bp_wall / 3600:.2f} h wall)", flush=True)

    # ---- ladder arms (independent streams) ----
    rung_arms: dict[str, dict[str, object]] = {}
    rung_verdicts: dict[str, dict[str, object]] = {}
    for rung in RUNGS:
        print(f"{rung}: decoding own stream ...", flush=True)
        meta = run_beam(
            cdir, rung, f"shots_{rung}.bin", f"{rung}_predictions.bin", rung,
            threads, str(arms[rung]["sha256"]))  # type: ignore[index]
        failed, per_obs_r = failure_mask_from_predictions(
            cdir / f"shots_{rung}.bin", cdir / f"{rung}_predictions.bin")
        mask_sha = write_failmask(cdir / f"{rung}_failmask.bin", failed)
        arm = arm_summary(rung, shots, int(failed.sum()), per_obs_r.astype(int).tolist(),
                          meta, f"{rung}_failmask.bin", mask_sha)
        arm["source_shots_sha256"] = str(arms[rung]["sha256"])  # type: ignore[index]
        verdict = ratio_verdict(arm, bposd, RUNGS[rung]["band"],
                                derive_seed(BASE_SEED, "gateB-gb5a-ratio-bootstrap", rung, shots))
        campaign.append_result(arm)
        rung_arms[rung] = arm
        rung_verdicts[rung] = verdict
        print(f"  {rung}: {arm['failures']}/{shots} failures -> {verdict['outcome']} "
              f"(ratio {verdict['beam_over_bposd_ratio']}, CI {verdict['ratio_ci95']})", flush=True)

    # ---- paired instrument on the shared BP+OSD stream ----
    paired_masks: dict[str, np.ndarray] = {"bp30+osd": bp_mask}
    paired_meta: dict[str, object] = {}
    for rung in RUNGS:
        print(f"{rung}: decoding shared paired stream ...", flush=True)
        meta = run_beam(
            cdir, rung, "shots_bposd.bin", f"paired_{rung}_predictions.bin",
            f"paired_{rung}", threads, bp_sha)
        failed, _po = failure_mask_from_predictions(
            cdir / "shots_bposd.bin", cdir / f"paired_{rung}_predictions.bin")
        write_failmask(cdir / f"paired_{rung}_failmask.bin", failed)
        paired_masks[rung] = failed
        paired_meta[rung] = meta
        print(f"  paired {rung}: {int(failed.sum())}/{shots} failures", flush=True)

    names = ["bp30+osd", *RUNGS]
    discordance = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            fa, fb = paired_masks[a], paired_masks[b]
            n_ab = int((fa & ~fb).sum())
            n_ba = int((fb & ~fa).sum())
            discordance[f"{a}|{b}"] = {
                "a_fails_b_ok": n_ab, "b_fails_a_ok": n_ba,
                "both_fail": int((fa & fb).sum()),
                "mcnemar_exact_p": mcnemar_exact(n_ab, n_ba),
            }
    identity = {}
    for i, a in enumerate(list(RUNGS)):
        for b in list(RUNGS)[i + 1:]:
            sa = sha256_file(cdir / f"paired_{a}_predictions.bin")
            sb = sha256_file(cdir / f"paired_{b}_predictions.bin")
            identity[f"{a}|{b}"] = {"identical_predictions": sa == sb,
                                    "sha256_a": sa, "sha256_b": sb}
    paired = {
        "shared_stream": "shots_bposd.bin",
        "shared_stream_sha256": bp_sha,
        "shots": shots,
        "failures": {k: int(v.sum()) for k, v in paired_masks.items()},
        "discordance": discordance,
        "prediction_identity": identity,
        "decode_meta": paired_meta,
        "interpretation": (
            "published ladder implies beam32 fixes ~82% and beam64 ~94% of beam8 "
            "failures; zero discordance is a sample identity refuting that, any "
            "nonzero discordance measures the in-harness width benefit"
        ),
    }

    summary = {
        "gate": "B-gb5a-ladder+paired",
        "p": P, "basis": BASIS, "shots_per_arm": shots,
        "bposd": bposd,
        "rung_arms": rung_arms,
        "rung_verdicts": rung_verdicts,
        "ladder_verdict": ladder_verdict(rung_verdicts),
        "paired_instrument": paired,
        "rate_semantics": {
            "raw_logical_failure_rate": "any-observable failure / shots",
            "per_round_ler": "raw_logical_failure_rate / 12",
            "ratio": "beam raw rate / BP+OSD raw rate; identical under per-round convention",
            "ratio_ci": "independent-binomial percentile; positive/zero draws retained as +inf",
        },
        "dependence_disclosure": (
            "all three rung ratios share one BP+OSD denominator sample, so the "
            "three band decisions are statistically dependent"
        ),
        "scope": "terminal Gate-B ladder decision per Revisions GB5/GB5a",
    }
    campaign.close(summary)
    print(json.dumps({"campaign": str(cdir),
                      "ladder": summary["ladder_verdict"],
                      "paired_failures": paired["failures"]}, indent=2))
    return cdir


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shots", type=int, default=20_000_000)
    ap.add_argument("--shards", type=int, default=240)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--threads", type=int, default=26)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    if args.smoke:
        if args.shots == 20_000_000:
            args.shots = 20_000
        if args.shards == 240:
            args.shards = 4
        if args.workers == 24:
            args.workers = 4
        if args.threads == 26:
            args.threads = 4
    try:
        os.nice(5)
    except OSError:
        pass
    run(args)


if __name__ == "__main__":
    main()
