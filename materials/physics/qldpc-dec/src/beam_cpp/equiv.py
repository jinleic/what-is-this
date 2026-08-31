"""beam8 equivalence harness: shared shots -> Python beam8 vs C++ beam8.

Generates seeded shot sets from the committed circuits, decodes them with BOTH
arms on identical shots, and diffs per-shot logical-observable predictions.

Usage:
  python equiv.py --p 3e-3 --shots 1000 --seed 20260830 [--cpp-only]
                  [--cpp-repeats N] [--out scratch/equiv_p3e-3.json]

Artifact caching (all under src/beam_cpp/scratch/):
  shots_<tag>.bin      shared bit-packed shots (stim u64-LE packing)
  py_preds_<tag>.npy   Python arm predictions (n x nobs int64)
  cpp_preds_<tag>.bin  C++ arm raw output (see README for format)
"""
from __future__ import annotations

import argparse
import json
import struct
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import stim

HERE = Path(__file__).resolve().parent
TARGET = HERE.parent.parent.parent            # qldpc-dec/
SRC = TARGET / "src"
sys.path.insert(0, str(SRC))

from qldpc_dec.beam_search import BeamSearchBatchDecoder  # noqa: E402
from qldpc_dec.circuits import circuit_sha, dem_sha, load_circuit  # noqa: E402

STIM_VERSION = stim.__version__


# ---------------------------------------------------------------------------
# shared shot generation (bit-packed via struct: stim packs u64 little-endian)
# ---------------------------------------------------------------------------

def pack_bits(det: np.ndarray) -> np.ndarray:
    """(shots, nbits) bool -> (shots, ceil(nbits/64)*8) uint8 in stim's u64-LE
    packing (word w little-endian, bit b of word = column 64*w+b)."""
    shots, nbits = det.shape
    n64 = (nbits + 63) // 64
    padded = np.zeros((shots, n64 * 64), dtype=bool)
    padded[:, :nbits] = det
    words = padded.reshape(shots, n64, 64)
    weights = (np.uint64(1) << np.arange(64, dtype=np.uint64)).reshape(1, 1, 64)
    u64 = (words.astype(np.uint64) * weights).sum(axis=2, dtype=np.uint64)
    return np.ascontiguousarray(u64.astype("<u8")).view(np.uint8).reshape(shots, n64 * 8)


def sample_shots(circuit: stim.Circuit, shots: int, seed: int):
    """Return (packed_det, packed_obs) uint8 arrays + (ndet, nobs)."""
    sampler = circuit.compile_detector_sampler(seed=seed)
    det, obs = sampler.sample(shots, separate_observables=True)
    det = np.ascontiguousarray(det, dtype=bool)
    obs = np.ascontiguousarray(obs, dtype=bool)
    return pack_bits(det), pack_bits(obs), det.shape[1], obs.shape[1]


def save_shots(path: Path, packed_det: np.ndarray, packed_obs: np.ndarray,
               ndet: int, nobs: int) -> None:
    with open(path, "wb") as f:
        f.write(struct.pack("<QQQ", ndet, nobs, packed_det.shape[0]))
        f.write(np.ascontiguousarray(packed_det, dtype=np.uint8).tobytes())
        f.write(np.ascontiguousarray(packed_obs, dtype=np.uint8).tobytes())


# ---------------------------------------------------------------------------
# arms
# ---------------------------------------------------------------------------

def unpack_bits(packed: np.ndarray, nbits: int) -> np.ndarray:
    """Inverse of pack_bits: (shots, nbytes) uint8 -> (shots, nbits) bool."""
    shots = packed.shape[0]
    n64 = (nbits + 63) // 64
    u64 = packed.reshape(shots, n64 * 8).view("<u8") if False else \
        np.ascontiguousarray(packed).reshape(shots, n64, 8).view("<u8").reshape(shots, n64)
    bits = ((u64[:, :, None] >> np.arange(64, dtype=np.uint64)[None, None, :]) &
            np.uint64(1)).astype(bool)
    return bits.reshape(shots, n64 * 64)[:, :nbits]


def run_python_arm(dem, shots_bin: Path, ndet: int, tag: str) -> dict:
    from qldpc_dec.beam_search import BeamSearchBatchDecoder
    raw = shots_bin.read_bytes()
    hdr_ndet, hdr_nobs, n = struct.unpack("<QQQ", raw[:24])
    assert hdr_ndet == ndet
    packed = np.frombuffer(raw, dtype=np.uint8, offset=24)
    det_bytes_row = (hdr_ndet + 63) // 64 * 8   # u64-word padded rows
    obs_bytes_row = (hdr_nobs + 63) // 64 * 8
    assert len(packed) == n * (det_bytes_row + obs_bytes_row)
    det_packed = packed[: n * det_bytes_row].reshape(n, det_bytes_row)
    det = unpack_bits(np.ascontiguousarray(det_packed), hdr_ndet)
    dec = BeamSearchBatchDecoder(dem, beam_width=8, initial_iters=30,
                                 iters_per_round=20, max_rounds=10, num_results=1)
    t0 = time.perf_counter()
    preds = dec.decode_batch(det)
    wall = time.perf_counter() - t0
    np.save(HERE / "scratch" / f"py_preds_{tag}.npy", preds.astype(np.int64))
    return {"ms_per_shot": 1e3 * wall / n, "shots": n}


def run_cpp_arm(binpath: Path, dem_path: Path, lam_path: Path, shots_path: Path,
                tag: str, nobs: int, nshots: int, repeats: int = 1,
                threads: int = 1) -> dict:
    out_bin = HERE / "scratch" / f"cpp_preds_{tag}.bin"
    cmd = [str(binpath), "--dem", str(dem_path), "--lam", str(lam_path),
           "--shots", str(shots_path), "--out", str(out_bin)]
    if repeats > 1:
        cmd.append(f"--repeats={repeats}")
    if threads > 1:
        cmd.append(f"--threads={threads}")
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if r.returncode != 0:
        raise RuntimeError(f"C++ decoder failed:\n{r.stdout}\n{r.stderr}")
    sys.stderr.write(r.stderr)
    raw = out_bin.read_bytes()
    magic, hdr_nobs, hdr_n = struct.unpack("<III", raw[:12])
    assert (magic, hdr_nobs, hdr_n) == (0x31503842, nobs, nshots), \
        f"C++ output header mismatch: {magic:#x} {hdr_nobs} {hdr_n} vs nobs={nobs} n={nshots}"
    words = np.frombuffer(raw, dtype="<u8", offset=12).reshape(nshots, (nobs + 63) // 64)
    bits = ((words[:, :, None] >> np.arange(64, dtype=np.uint64)[None, None, :]) &
            np.uint64(1)).astype(np.int64)
    preds = bits.reshape(nshots, -1)[:, :nobs]
    np.save(HERE / "scratch" / f"cpp_preds_{tag}.npy", preds)
    # decode time: subprocess wall includes DEM parse; the binary prints its
    # own decode-only wall to stderr -- parse that for the timing number.
    dec_ms = None
    for line in r.stderr.splitlines():
        if "ms/shot" in line:
            dec_ms = float(line.split("=>")[1].replace("ms/shot", "").strip())
    return {"ms_per_shot_cpp": dec_ms, "subprocess_s": wall}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=float, default=3e-3)
    ap.add_argument("--shots", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260830)
    ap.add_argument("--cpp-only", action="store_true",
                    help="reuse cached Python predictions; decode with C++ only")
    ap.add_argument("--cpp-repeats", type=int, default=1,
                    help="repeat C++ decode this many times to cut timing noise")
    ap.add_argument("--cpp-threads", type=int, default=1)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.p == 3e-3:
        # no upstream circuit at p=3e-3; the companion's committed DERIVED
        # rescale artifact (sha ea2de75c...) is the exact circuit the frozen
        # beam8 2362 ms/shot reference used -- load it directly.
        derived = HERE / ".." / ".." / "circuits" / \
            "BB_144_144_12_memory_Z_p0.003_sr12_derived_p1e-3_rescale.stim"
        circuit = stim.Circuit.from_file(derived)
    else:
        circuit = load_circuit(args.p, "Z")
    dem = circuit.detector_error_model(decompose_errors=True,
                                       ignore_decomposition_failures=True)
    ndet, nobs = dem.num_detectors, dem.num_observables
    tag = f"p{args.p:g}_n{args.shots}_s{args.seed}"
    scratch = HERE / "scratch"
    scratch.mkdir(exist_ok=True)

    # --- common artifacts ---
    dem_path = scratch / f"dem_{tag}.dem"
    dem_path.write_text(str(dem))
    lam_path = scratch / f"lam_{tag}.npy"
    if not lam_path.exists():
        from qldpc_dec.dem_matrices import dem_to_matrices
        _H, _A, lam = dem_to_matrices(dem, merge=False)
        np.save(lam_path, lam)

    shots_path = scratch / f"shots_{tag}.bin"
    if not shots_path.exists():
        packed_det, packed_obs, ndet_s, nobs_s = sample_shots(circuit, args.shots,
                                                              args.seed)
        save_shots(shots_path, packed_det, packed_obs, ndet_s, nobs_s)

    py_cache = scratch / f"py_preds_{tag}.npy"
    py = {}
    if not args.cpp_only:
        py = run_python_arm(dem, shots_path, ndet, tag)
    else:
        assert py_cache.exists(), "--cpp-only needs cached py_preds"
        py = {"ms_per_shot": float("nan"), "shots": args.shots, "note": "cached"}

    cpp = run_cpp_arm(HERE / "beam8_cpp", dem_path, lam_path, shots_path, tag,
                      nobs, args.shots, repeats=max(1, args.cpp_repeats),
                      threads=args.cpp_threads)

    py_preds = np.load(py_cache)
    cpp_preds = np.load(scratch / f"cpp_preds_{tag}.npy")
    assert py_preds.shape == cpp_preds.shape, (py_preds.shape, cpp_preds.shape)
    mism = np.where((py_preds != cpp_preds).any(axis=1))[0]
    speedup = (py["ms_per_shot"] / cpp["ms_per_shot_cpp"]
               if py.get("ms_per_shot") == py.get("ms_per_shot") else None)
    sha = None
    if args.p == 3e-3:
        import hashlib
        derived = (HERE / ".." / ".." / "circuits" /
                   "BB_144_144_12_memory_Z_p0.003_sr12_derived_p1e-3_rescale.stim").resolve()
        sha = hashlib.sha256(derived.read_bytes()).hexdigest()
    else:
        sha = circuit_sha(args.p, "Z")
    report = {
        "p": args.p, "shots": args.shots, "seed": args.seed,
        "ndet": ndet, "nobs": nobs,
        "circuit_sha256": sha,
        "dem_sha256": dem_sha(args.p, "Z") if args.p != 3e-3 else "(derived circuit; compute via circuits.py after DEM committed)",
        "python_ms_per_shot": py.get("ms_per_shot"),
        "cpp_ms_per_shot": cpp["ms_per_shot_cpp"],
        "cpp_speedup_vs_python": speedup,
        "cpp_threads": args.cpp_threads,
        "python_failures": int((py_preds != 0).any(axis=1).sum()) if py_preds.shape[1] else None,
        "mismatched_shots": int(len(mism)),
        "mismatched_shot_ids": mism[:100].tolist(),
    }
    out = args.out or str(scratch / f"equiv_{tag}.json")
    Path(out).write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))
    if len(mism):
        print(f"[DIVERGENCE] {len(mism)}/{args.shots} shots differ; first ids: "
              f"{mism[:10].tolist()}")
    else:
        print(f"[MATCH] {args.shots}/{args.shots} shots identical.")


if __name__ == "__main__":
    main()
