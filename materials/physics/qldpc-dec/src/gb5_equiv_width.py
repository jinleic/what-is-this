"""Bit-exactness check of the C++ port at ARBITRARY Table-I beam parameters.

GB1 established 0/1000 mismatches for the beam8 parameter set only, and the
C++ port silently ignored its own --beam-width/--initial-iters/
--iters-per-round/--max-rounds flags until the 2026-08-31 wiring fix. This
runner re-establishes the equivalence contract per rung:

    python  BeamSearchBatchDecoder(dem, **params)
    cpp     <binary> --beam-width=.. --initial-iters=.. --iters-per-round=..
                     --max-rounds=.. --num-results=..

on one shared seeded shot stream, and writes a frozen JSON report with the
binary sha256, the parameters, the mismatch count, and the shot ids of any
divergence. Zero mismatches is the prerequisite for using the C++ arm.

Revision GB9 (2026-09-02): rung ``beam64_32res`` (Table 1 beam64_32res_640iters,
num_results=32) and ``--binary`` (default ``beam8_cpp`` = the frozen GB5a/GB6/GB7
binary; ``beam_nr_cpp`` carries --num-results). Reports for a non-default binary
get the binary name in their filename so the six frozen reports are never
overwritten. ``--py-workers`` shards the (slow, pure-numpy) reference decode
across processes; per-shot results are independent so sharding is exact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import struct
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import stim

TARGET = Path(__file__).resolve().parents[1]
SRC = TARGET / "src"
BEAM_DIR = SRC / "beam_cpp"
sys.path.insert(0, str(SRC))

from qldpc_dec.beam_search import BeamSearchBatchDecoder  # noqa: E402
from qldpc_dec.circuits import circuit_sha, dem_sha, load_circuit, load_dem  # noqa: E402
from qldpc_dec.dem_matrices import dem_to_matrices  # noqa: E402
from qldpc_dec.seeds import derive_seed  # noqa: E402

RUNG_PARAMS = {
    "beam8": dict(beam_width=8, initial_iters=30, iters_per_round=20, max_rounds=10,
                  num_results=1),
    "beam32": dict(beam_width=32, initial_iters=40, iters_per_round=30, max_rounds=10,
                   num_results=1),
    "beam64": dict(beam_width=64, initial_iters=40, iters_per_round=30, max_rounds=20,
                   num_results=1),
    # arXiv:2512.07057 Table 1 row 4: the configuration behind the 17x claim.
    "beam64_32res": dict(beam_width=64, initial_iters=40, iters_per_round=30,
                         max_rounds=20, num_results=32),
}
DEFAULT_BINARY = "beam8_cpp"
OUT_DIR = BEAM_DIR / "evidence"

_WORKER_DECODER: BeamSearchBatchDecoder | None = None


def _worker_init(p: float, params: dict) -> None:
    # Rebuild the DEM through the same deterministic loader as the parent
    # (no text round-trip), so every worker's lam is bit-identical.
    global _WORKER_DECODER
    _WORKER_DECODER = BeamSearchBatchDecoder(load_dem(p, "Z"), **params)


def _worker_decode(chunk: np.ndarray) -> np.ndarray:
    assert _WORKER_DECODER is not None
    return np.asarray(_WORKER_DECODER.decode_batch(chunk), dtype=np.int64)


def python_predictions(p: float, dem: stim.DetectorErrorModel, params: dict,
                       det: np.ndarray, workers: int) -> np.ndarray:
    if workers <= 1:
        return np.asarray(BeamSearchBatchDecoder(dem, **params).decode_batch(det),
                          dtype=np.int64)
    import concurrent.futures
    import multiprocessing
    chunks = [c for c in np.array_split(np.arange(len(det)), min(workers * 4, len(det)))
              if len(c)]
    ctx = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers, mp_context=ctx,
            initializer=_worker_init, initargs=(p, params)) as pool:
        parts = list(pool.map(_worker_decode, [det[c] for c in chunks]))
    return np.concatenate(parts, axis=0)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_shots(path: Path, circuit: stim.Circuit, shots: int, seed: int) -> tuple[int, int]:
    det_p, obs_p = circuit.compile_detector_sampler(seed=seed).sample(
        shots, separate_observables=True, bit_packed=True
    )
    ndet, nobs = circuit.num_detectors, circuit.num_observables
    db = (ndet + 63) // 64 * 8
    ob = (nobs + 63) // 64 * 8
    det = np.zeros((shots, db), dtype=np.uint8)
    det[:, : det_p.shape[1]] = det_p
    obs = np.zeros((shots, ob), dtype=np.uint8)
    obs[:, : obs_p.shape[1]] = obs_p
    with path.open("wb") as f:
        f.write(struct.pack("<QQQ", ndet, nobs, shots))
        f.write(det.tobytes())
        f.write(obs.tobytes())
    return ndet, nobs


def load_det(path: Path) -> np.ndarray:
    raw = path.read_bytes()
    ndet, nobs, shots = struct.unpack("<QQQ", raw[:24])
    db = (ndet + 63) // 64 * 8
    det_p = np.frombuffer(raw, dtype=np.uint8, offset=24, count=shots * db).reshape(shots, db)
    return np.unpackbits(det_p, axis=1, bitorder="little")[:, :ndet]


def load_preds(path: Path, nobs: int, shots: int) -> np.ndarray:
    raw = path.read_bytes()
    magic, hdr_nobs, hdr_shots = struct.unpack("<III", raw[:12])
    assert (magic, hdr_nobs, hdr_shots) == (0x31503842, nobs, shots), (magic, hdr_nobs, hdr_shots)
    words = np.frombuffer(raw, dtype="<u8", offset=12).reshape(shots, (nobs + 63) // 64)
    bits = (
        (words[:, :, None] >> np.arange(64, dtype=np.uint64)[None, None, :]) & np.uint64(1)
    ).astype(np.int64)
    return bits.reshape(shots, -1)[:, :nobs]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rung", required=True, choices=list(RUNG_PARAMS))
    ap.add_argument("--p", type=float, required=True)
    ap.add_argument("--shots", type=int, required=True)
    ap.add_argument("--binary", default=DEFAULT_BINARY,
                    help="binary name under src/beam_cpp (default: frozen beam8_cpp)")
    ap.add_argument("--cpp-threads", type=int, default=1)
    ap.add_argument("--py-workers", type=int, default=1)
    args = ap.parse_args()
    beam_bin = BEAM_DIR / args.binary
    if not beam_bin.is_file():
        raise FileNotFoundError(beam_bin)
    if RUNG_PARAMS[args.rung]["num_results"] > 1 and args.binary == DEFAULT_BINARY:
        raise SystemExit("the frozen beam8_cpp has no --num-results; use --binary beam_nr_cpp")
    report_suffix = "" if args.binary == DEFAULT_BINARY else f"_{args.binary}"

    params = RUNG_PARAMS[args.rung]
    circuit = load_circuit(args.p, "Z")
    dem = load_dem(args.p, "Z")
    OUT_DIR.mkdir(exist_ok=True)
    tag = f"{args.rung}_p{args.p:g}_n{args.shots}"
    dem_path = OUT_DIR / f"gb5eq_{tag}.dem"
    dem_path.write_text(str(dem))
    _H, _A, lam = dem_to_matrices(dem, merge=False)
    lam_path = OUT_DIR / f"gb5eq_{tag}_lam.npy"
    np.save(lam_path, lam)

    seed = derive_seed(20260829, "gb5-equiv", args.rung, args.p, args.shots)
    shots_path = OUT_DIR / f"gb5eq_{tag}_shots.bin"
    ndet, nobs = write_shots(shots_path, circuit, args.shots, seed)

    det = load_det(shots_path)
    cache_path = OUT_DIR / f"gb5eq_{tag}_py.npy"
    if cache_path.exists():
        py_preds = np.load(cache_path)
        if py_preds.shape[0] != args.shots:
            raise ValueError(f"cached python predictions have wrong shape: {py_preds.shape}")
        py_ms = None  # cached from an earlier identical (rung, p, shots, seed) run
    else:
        t0 = time.perf_counter()
        py_preds = python_predictions(args.p, dem, params, det, args.py_workers)
        py_ms = 1000.0 * (time.perf_counter() - t0) / args.shots
        np.save(cache_path, py_preds)

    out_path = OUT_DIR / f"gb5eq_{tag}{report_suffix}_cpp.bin"
    cmd = [
        str(beam_bin),
        "--dem", str(dem_path),
        "--lam", str(lam_path),
        "--shots", str(shots_path),
        "--out", str(out_path),
        f"--beam-width={params['beam_width']}",
        f"--initial-iters={params['initial_iters']}",
        f"--iters-per-round={params['iters_per_round']}",
        f"--max-rounds={params['max_rounds']}",
        "--threads", str(args.cpp_threads),
    ]
    if params["num_results"] > 1 or args.binary != DEFAULT_BINARY:
        cmd.append(f"--num-results={params['num_results']}")
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    cpp_wall = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr)
    cfg_line = next(l for l in proc.stderr.splitlines() if l.startswith("config"))
    cpp_ms = None
    for line in proc.stderr.splitlines():
        if "ms/shot" in line:
            cpp_ms = float(line.split("=>")[1].replace("ms/shot", "").strip())

    cpp_preds = load_preds(out_path, nobs, args.shots)
    mism = np.where((np.asarray(py_preds, dtype=np.int64) != cpp_preds).any(axis=1))[0]
    report = {
        "rung": args.rung,
        "parameters": dict(params),
        "binary": args.binary,
        "cpp_threads": args.cpp_threads,
        "python_workers": args.py_workers,
        "cpp_config_line": cfg_line.strip(),
        "p": args.p,
        "basis": "Z",
        "shots": args.shots,
        "seed": seed,
        "ndet": int(ndet),
        "nobs": int(nobs),
        "circuit_sha256": circuit_sha(args.p, "Z"),
        "dem_sha256": dem_sha(args.p, "Z"),
        "beam_binary_sha256": sha256_file(beam_bin),
        "beam_source_sha256": sha256_file(BEAM_DIR / "beam8.cpp"),
        "python_ms_per_shot": py_ms,
        "cpp_ms_per_shot": cpp_ms,
        "cpp_subprocess_s": cpp_wall,
        "mismatched_shots": int(len(mism)),
        "mismatched_shot_ids": mism[:100].tolist(),
        "shots_sha256": sha256_file(shots_path),
        "cpp_predictions_sha256": sha256_file(out_path),
    }
    dest = OUT_DIR / f"gb5eq_{tag}{report_suffix}.json"
    dest.write_text(json.dumps(report, indent=2))
    verdict = "MATCH" if not len(mism) else "DIVERGENCE"
    print(f"[{verdict}] {tag}: {args.shots - len(mism)}/{args.shots} identical; report {dest.name}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
