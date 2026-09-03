"""GB9 - the paper's actual 17x claim: beam64_32res_640iters (num_results=32).

pre_statement.md Revision GB9 (recorded before any sample). Decode beam8 and
beam64_32res (64, 40, 30, 20, num_results=32; arXiv:2512.07057 Table 1 row 4,
Section III "17x") on the SAME 1e8-shot bp30+osd stream that GB7 decoded
(sampling_seed(20260829, 1e-3, "Z", "bp30+osd")). Hard voids:

  * the resampled stream must equal the frozen GB7 shots_bposd.bin byte-for-byte;
  * the beam8 predictions (new binary beam_nr_cpp, --num-results=1) must equal
    the frozen GB7 beam8_predictions.bin byte-for-byte (1e8-shot K=1 regression);
  * the beam8 failure mask must equal the frozen GB7 mask.

Decoding is shard-wise (SHARD shots per ordinary shots file, same binary,
predictions concatenated): per-shot results are independent, so the assembled
file is byte-identical to a single run, and the run is resumable through the
idempotent per-shard meta files. Decision rule and verdict mapping: see
pre_statement.md Revision GB9 and `decide()` below.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
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
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gateb_gb5a_ladder import (  # noqa: E402
    BASE_SEED, BASIS, BEAM_DIR, BEAM_SEARCH_SOURCE, BEAM_SOURCE, BOOTSTRAP_SOURCE,
    CHUNK, DEM_MATRICES_SOURCE, RUNGS, check_predictions_header,
    failure_mask_from_predictions, mcnemar_exact, read_failmask, read_shots_header,
    sha256_file, write_failmask, write_shots,
)
from gateb_gb7_paired import (  # noqa: E402
    frozen_checksums, paired_cells, verify_frozen_file,
)
from qldpc_dec.bootstrap import wilson_interval  # noqa: E402
from qldpc_dec.circuits import circuit_sha, dem_sha, load_circuit, load_dem  # noqa: E402
from qldpc_dec.dem_matrices import dem_to_matrices  # noqa: E402
from qldpc_dec.runner import Campaign  # noqa: E402
from qldpc_dec.seeds import derive_seed, sampling_seed  # noqa: E402

P = 1e-3
ROUNDS = 12
SHOTS = 100_000_000
SHARD = 5_000_000
BOOTSTRAP_DRAWS = 100_000
MCNEMAR_ALPHA = 0.05
PUBLISHED_FACTOR = 17.0
BAND_FACTOR = (PUBLISHED_FACTOR / 1.3, PUBLISHED_FACTOR * 1.3)
GB7_CAMPAIGN_ID = "20260901T145247Z_b7ea9ac4_ac3f6689e03b"
GB7_CAMPAIGN = TARGET / "campaigns" / GB7_CAMPAIGN_ID
BEAM_NR_BINARY = BEAM_DIR / "beam_nr_cpp"
EVIDENCE = BEAM_DIR / "evidence"
SIZING = SRC / "evidence" / "gb9_sizing.json"
ARMS = MappingProxyType({
    "beam8": MappingProxyType({
        **{k: RUNGS["beam8"][k] for k in (
            "beam_width", "initial_iters", "iters_per_round", "max_rounds")},
        "num_results": 1,
        "paper_name": "beam8_230iters",
    }),
    "beam64_32res": MappingProxyType({
        "beam_width": 64, "initial_iters": 40, "iters_per_round": 30, "max_rounds": 20,
        "num_results": 32,
        "paper_name": "beam64_32res_640iters",
    }),
})
# Equivalence gate points for beam64_32res (binary beam_nr_cpp) and the six
# K=1 regression points whose C++ predictions must be byte-identical to the
# frozen beam8_cpp reports.
GATE_32RES = ((0.001, 2000), (0.003, 300))
GATE_K1 = (("beam8", 0.001, 2000), ("beam32", 0.001, 2000), ("beam64", 0.001, 2000),
           ("beam8", 0.003, 300), ("beam32", 0.003, 300), ("beam64", 0.003, 300))


# ----------------------------------------------------------------- launch gates

def equivalence_evidence_nr() -> dict[str, object]:
    binary_sha = sha256_file(BEAM_NR_BINARY)
    out: dict[str, object] = {"binary_sha256": binary_sha, "beam64_32res": [], "k1_regression": []}
    for p, n in GATE_32RES:
        path = EVIDENCE / f"gb5eq_beam64_32res_p{p:g}_n{n}_beam_nr_cpp.json"
        if not path.is_file():
            raise FileNotFoundError(f"missing 32res equivalence report {path}")
        rep = json.loads(path.read_text())
        if rep["mismatched_shots"] != 0:
            raise AssertionError(f"32res equivalence FAILED: {path}")
        if rep["beam_binary_sha256"] != binary_sha:
            raise AssertionError(f"{path} is for a different binary")
        want = {k: ARMS["beam64_32res"][k] for k in (
            "beam_width", "initial_iters", "iters_per_round", "max_rounds", "num_results")}
        if {k: rep["parameters"].get(k) for k in want} != want:
            raise AssertionError(f"{path} parameters {rep['parameters']} != {want}")
        out["beam64_32res"].append({
            "path": str(path.relative_to(TARGET)), "p": p, "shots": n,
            "mismatched_shots": 0, "sha256": sha256_file(path)})
    for rung, p, n in GATE_K1:
        new = EVIDENCE / f"gb5eq_{rung}_p{p:g}_n{n}_beam_nr_cpp.json"
        old = EVIDENCE / f"gb5eq_{rung}_p{p:g}_n{n}.json"
        rn, ro = json.loads(new.read_text()), json.loads(old.read_text())
        if rn["mismatched_shots"] != 0 or ro["mismatched_shots"] != 0:
            raise AssertionError(f"K=1 equivalence FAILED: {new}")
        if rn["beam_binary_sha256"] != binary_sha:
            raise AssertionError(f"{new} is for a different binary")
        if rn["cpp_predictions_sha256"] != ro["cpp_predictions_sha256"]:
            raise AssertionError(f"K=1 regression: {rung} p={p} predictions differ old vs new binary")
        out["k1_regression"].append({
            "rung": rung, "p": p, "shots": n,
            "cpp_predictions_sha256": rn["cpp_predictions_sha256"],
            "identical_to_frozen_binary": True})
    return out


def validate_gb7_campaign() -> dict[str, object]:
    status = json.loads((GB7_CAMPAIGN / "status.json").read_text())
    verdict = status.get("verdict")
    if not isinstance(verdict, str) or not verdict.startswith("FROZEN-"):
        raise ValueError(f"GB7 campaign is not evidentially frozen: {status}")
    checksums = frozen_checksums(GB7_CAMPAIGN)
    required = ["manifest.json", "summary.json", "shots_bposd.bin",
                "beam8_predictions.bin", "beam64_predictions.bin",
                "beam8_failmask.bin", "beam64_failmask.bin"]
    verified = {name: verify_frozen_file(GB7_CAMPAIGN, name, checksums) for name in required}
    summary = json.loads((GB7_CAMPAIGN / "summary.json").read_text())
    if int(summary.get("shots_shared", -1)) != SHOTS:
        raise ValueError("GB7 campaign does not contain exactly 1e8 shared shots")
    return {
        "run_id": GB7_CAMPAIGN_ID,
        "terminal_verdict": verdict,
        "sha256s_txt_sha256": sha256_file(GB7_CAMPAIGN / "sha256s.txt"),
        "verified_files": verified,
        "beam8_failures": summary["arms"]["beam8"]["failures"],
        "beam64_failures": summary["arms"]["beam64"]["failures"],
        "beam64_decision": summary["decisions"]["beam64"],
    }


# ------------------------------------------------------------------ sharding

def write_shards(cdir: Path, shots_file: str, nshards: int) -> list[dict[str, object]]:
    """Slice the shared shots file into ordinary per-shard shots files."""
    src = cdir / shots_file
    ndet, nobs, shots, det_bytes, obs_bytes = read_shots_header(src)
    if shots % nshards:
        raise ValueError(f"{shots} shots cannot be split into {nshards} equal shards")
    per = shots // nshards
    sdir = cdir / "shards"
    sdir.mkdir(exist_ok=True)
    inv_path = sdir / "inventory.json"
    if inv_path.is_file():
        inv = json.loads(inv_path.read_text())
        if inv.get("source_sha256") == sha256_file(src) and len(inv.get("shards", [])) == nshards:
            for s in inv["shards"]:
                if (sdir / s["file"]).stat().st_size != s["bytes"]:
                    raise ValueError(f"shard {s['file']} size drift")
            return inv["shards"]
        raise ValueError("stale shard inventory")
    shards = []
    with src.open("rb") as f:
        for k in range(nshards):
            name = f"shots_{k:02d}.bin"
            path = sdir / name
            with path.open("wb") as g:
                g.write(struct.pack("<QQQ", ndet, nobs, per))
                f.seek(24 + k * per * det_bytes)
                remaining = per * det_bytes
                while remaining:
                    block = f.read(min(1 << 24, remaining))
                    g.write(block)
                    remaining -= len(block)
                f.seek(24 + shots * det_bytes + k * per * obs_bytes)
                remaining = per * obs_bytes
                while remaining:
                    block = f.read(min(1 << 24, remaining))
                    g.write(block)
                    remaining -= len(block)
                g.flush()
                os.fsync(g.fileno())
            shards.append({"index": k, "file": name, "start": k * per, "shots": per,
                           "bytes": path.stat().st_size, "sha256": sha256_file(path)})
            print(f"shard {k:02d}/{nshards}: {name} written", flush=True)
    inv_path.write_text(json.dumps({"source": shots_file, "source_sha256": sha256_file(src),
                                    "shards": shards}, indent=2))
    return shards


def run_beam_nr(cdir: Path, arm: str, shots_path: Path, out_path: Path, log_tag: str,
                threads: int, source_shots_sha256: str) -> dict[str, object]:
    """Decode one shot file with one arm's pinned parameters on beam_nr_cpp (idempotent)."""
    params = ARMS[arm]
    expected_cfg = (
        f"config: beam_width={params['beam_width']} initial_iters={params['initial_iters']} "
        f"iters_per_round={params['iters_per_round']} max_rounds={params['max_rounds']} "
        f"num_results={params['num_results']}")
    binary_sha = sha256_file(BEAM_NR_BINARY)
    meta_path = out_path.with_name(out_path.name + ".meta.json")
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        valid = (
            out_path.is_file()
            and meta.get("predictions_sha256") == sha256_file(out_path)
            and meta.get("beam_binary_sha256") == binary_sha
            and meta.get("arm") == arm
            and meta.get("source_shots_sha256") == source_shots_sha256
            and meta.get("config_line") == expected_cfg
        )
        if valid:
            return meta
        raise ValueError(f"stale beam predictions or metadata for {out_path.name}")
    cmd = [
        str(BEAM_NR_BINARY),
        "--dem", str(cdir / "pinned.dem"),
        "--lam", str(cdir / "pinned_lam.npy"),
        "--shots", str(shots_path),
        "--out", str(out_path),
        f"--beam-width={params['beam_width']}",
        f"--initial-iters={params['initial_iters']}",
        f"--iters-per-round={params['iters_per_round']}",
        f"--max-rounds={params['max_rounds']}",
        f"--num-results={params['num_results']}",
        "--threads", str(threads),
    ]
    stderr_path = out_path.with_name(f"{log_tag}.stderr.log")
    t0 = time.perf_counter()
    with stderr_path.open("w") as se:
        rc = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=se)
    wall = time.perf_counter() - t0
    if rc != 0:
        raise RuntimeError(f"{arm} decode failed ({rc}):\n{stderr_path.read_text()[-4000:]}")
    stderr = stderr_path.read_text()
    cfg = next((l.strip() for l in stderr.splitlines() if l.startswith("config:")), "")
    if cfg != expected_cfg:
        raise AssertionError(f"binary ran the wrong configuration:\n  got {cfg}\n  want {expected_cfg}")
    stats = next((l.strip() for l in stderr.splitlines() if l.startswith("stats:")), "")
    ms_match = re.search(r"=>\s*([0-9.eE+-]+)\s*ms/shot", stderr)
    meta = {
        "arm": arm,
        "shots_file": shots_path.name,
        "predictions": out_path.name,
        "predictions_sha256": sha256_file(out_path),
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


def assemble_predictions(cdir: Path, arm: str, shards: list[dict[str, object]],
                         nobs: int, total: int) -> Path:
    """Concatenate shard predictions into one file with the full-run header."""
    out = cdir / f"{arm}_predictions.bin"
    row_bytes = ((nobs + 63) // 64) * 8
    with out.open("wb") as g:
        g.write(struct.pack("<III", 0x31503842, nobs, total))
        for s in shards:
            part = cdir / "shards" / f"{arm}_{s['index']:02d}_predictions.bin"
            check_predictions_header(part, nobs, int(s["shots"]))
            with part.open("rb") as f:
                f.seek(12)
                remaining = int(s["shots"]) * row_bytes
                while remaining:
                    block = f.read(min(1 << 24, remaining))
                    if not block:
                        raise ValueError(f"truncated shard predictions {part}")
                    g.write(block)
                    remaining -= len(block)
        g.flush()
        os.fsync(g.fileno())
    check_predictions_header(out, nobs, total)
    return out


def decode_arm(cdir: Path, arm: str, shards: list[dict[str, object]], threads: int,
               nobs: int, total: int = SHOTS) -> dict[str, object]:
    metas = []
    for s in shards:
        k = int(s["index"])
        shots_path = cdir / "shards" / str(s["file"])
        out_path = cdir / "shards" / f"{arm}_{k:02d}_predictions.bin"
        print(f"{arm}: shard {k:02d}/{len(shards)} ({s['shots']} shots) ...", flush=True)
        metas.append(run_beam_nr(cdir, arm, shots_path, out_path, f"{arm}_{k:02d}",
                                 threads, str(s["sha256"])))
    assembled = assemble_predictions(cdir, arm, shards, nobs, total)
    return {
        "shards": metas,
        "assembled": assembled.name,
        "predictions_sha256": sha256_file(assembled),
        "wall_s_total": float(sum(m["wall_s"] for m in metas)),
        "ms_per_shot_wall": 1000.0 * float(sum(m["wall_s"] for m in metas)) / total,
    }


# ------------------------------------------------------------------ decision

def paired_ratio_interval(cells: dict[str, int], shots: int) -> dict[str, object]:
    both = int(cells["both_fail"])
    ref_only = int(cells["reference_fails_candidate_ok"])
    candidate_only = int(cells["reference_ok_candidate_fails"])
    neither = int(cells["neither_fails"])
    counts = np.asarray([both, ref_only, candidate_only, neither], dtype=np.int64)
    if np.any(counts < 0) or int(counts.sum()) != shots:
        raise ValueError(f"invalid paired table: {counts.tolist()}")
    f8 = both + ref_only
    fr = both + candidate_only
    seed = derive_seed(BASE_SEED, "gb9-paired-ratio", "beam64_32res", shots)
    if f8 == 0:
        return {"method": "undefined_no_beam8_failures", "ratio": None, "ci95": [None, None],
                "beam8_failures": 0, "rung_failures": fr, "bootstrap_seed": seed}
    ratio = fr / f8
    if candidate_only == 0:
        lo, hi = wilson_interval(both, f8)
        method, draws, valid_draws = "Wilson on both_fail / beam8_failures (c=0)", None, None
    else:
        rng = np.random.default_rng(seed)
        samples = rng.multinomial(shots, counts / shots, size=BOOTSTRAP_DRAWS)
        den = samples[:, 0] + samples[:, 1]
        num = samples[:, 0] + samples[:, 2]
        valid = den > 0
        boot = num[valid] / den[valid]
        if boot.size == 0:
            raise RuntimeError("all paired bootstrap draws had zero beam8 failures")
        lo, hi = np.quantile(boot, [0.025, 0.975])
        method, draws, valid_draws = "four-cell paired multinomial percentile bootstrap", BOOTSTRAP_DRAWS, int(boot.size)
    return {"method": method, "ratio": float(ratio), "ci95": [float(lo), float(hi)],
            "beam8_failures": f8, "rung_failures": fr, "bootstrap_draws": draws,
            "valid_bootstrap_draws": valid_draws, "bootstrap_seed": seed}


def decide(cells8: dict[str, int], cells64: dict[str, int], shots: int,
           k1_ratio: float) -> dict[str, object]:
    """Pre-registered GB9 rule (pre_statement.md Revision GB9)."""
    b, c = int(cells8["reference_fails_candidate_ok"]), int(cells8["reference_ok_candidate_fails"])
    p_value = mcnemar_exact(b, c)
    width_effect = bool(p_value < MCNEMAR_ALPHA and b > c)
    interval = paired_ratio_interval(cells8, shots)
    lo, hi = interval["ci95"]
    ratio = interval["ratio"]
    reciprocal = 1.0 / PUBLISHED_FACTOR
    band_lo, band_hi = 1.0 / BAND_FACTOR[1], 1.0 / BAND_FACTOR[0]
    if ratio is None or lo is None or hi is None:
        excludes_one = published_excluded = False
        band = "BAND_UNDEFINED"
        factor, factor_ci, k1_excluded = None, [None, None], False
    else:
        excludes_one = bool(hi < 1.0 or lo > 1.0)
        published_excluded = bool(reciprocal < lo or reciprocal > hi)
        if band_lo <= lo and hi <= band_hi:
            band = "BAND_REPRODUCED"
        elif hi < band_lo or lo > band_hi:
            band = "BAND_NOT_REPRODUCED"
        else:
            band = "BAND_INCONCLUSIVE"
        factor = None if ratio == 0 else float(1.0 / ratio)
        factor_ci = [float(1.0 / hi), None if lo == 0 else float(1.0 / lo)]
        k1_excluded = bool(k1_ratio < lo or k1_ratio > hi)
    gap = bool(width_effect and excludes_one and published_excluded)
    satisfied = [label for label, ok in (
        ("WIDTH_EFFECT_CONFIRMED", width_effect),
        ("INTERVAL_EXCLUDES_ONE", excludes_one),
        ("PUBLISHED_FACTOR_EXCLUDED", published_excluded)) if ok]
    if gap:
        outcome = "GAP_CONFIRMED"
    elif satisfied:
        outcome = "NO_GAP_" + "_AND_".join(satisfied)
    else:
        outcome = "UNRESOLVED_AT_N"
    # Same-shot comparison against the frozen GB7 beam64_640iters (num_results=1).
    b64 = int(cells64["reference_fails_candidate_ok"])
    c64 = int(cells64["reference_ok_candidate_fails"])
    p64 = mcnemar_exact(b64, c64)
    improvement = bool(p64 < MCNEMAR_ALPHA and b64 > c64)
    if gap:
        verdict = "FROZEN-NEGATIVE"
    elif not published_excluded and improvement:
        verdict = "FROZEN-CERTIFIED"
    else:
        verdict = "FROZEN-INCONCLUSIVE"
    return {
        "rung": "beam64_32res",
        "paired_cells_vs_beam8": dict(cells8),
        "mcnemar_vs_beam8": {"exact_two_sided_p": p_value, "alpha": MCNEMAR_ALPHA,
                             "b_beam8_fails_rung_ok": b, "c_beam8_ok_rung_fails": c,
                             "width_effect_confirmed": width_effect},
        "paired_ratio_rung_over_beam8": interval,
        "interval_excludes_one": excludes_one,
        "published_factor": PUBLISHED_FACTOR,
        "published_reciprocal": reciprocal,
        "published_factor_excluded": published_excluded,
        "measured_factor_beam8_over_rung": factor,
        "measured_factor_ci95": factor_ci,
        "band_factor": list(BAND_FACTOR),
        "band_outcome": band,
        "gap_confirmed": gap,
        "outcome": outcome,
        "vs_beam64_640iters": {
            "paired_cells": dict(cells64),
            "frozen_gb7_beam64_ratio_over_beam8": k1_ratio,
            "k1_ratio_excluded_by_interval": k1_excluded,
            "mcnemar_exact_two_sided_p": p64,
            "b_beam64_fails_rung_ok": b64, "c_beam64_ok_rung_fails": c64,
            "improvement_over_k1_confirmed": improvement,
        },
        "terminal_verdict_by_prereg_mapping": verdict,
    }


# ----------------------------------------------------------------------- run

def run(args: argparse.Namespace) -> Path:
    run_dir = Path(args.run_dir).resolve()
    threads = int(args.threads)
    nshards = int(args.shards)
    if threads <= 0 or nshards <= 0 or SHOTS % nshards or (SHOTS // nshards) % CHUNK:
        raise ValueError("threads/shards must be positive; shards must be CHUNK-aligned")
    gb7 = validate_gb7_campaign()
    equivalence = equivalence_evidence_nr()
    seed = sampling_seed(BASE_SEED, P, BASIS, "bp30+osd")
    config = {
        "run": "gateB-gb9-beam64-32res",
        "p": P, "basis": BASIS, "rounds": ROUNDS, "base_seed": BASE_SEED,
        "sampling_seed": seed, "shots_shared": SHOTS, "shard_shots": SHOTS // nshards,
        "prefix_campaign_id": GB7_CAMPAIGN_ID, "prefix_shots": SHOTS,
        "arms": {arm: {k: v for k, v in ARMS[arm].items()} for arm in ARMS},
        "published_factor": PUBLISHED_FACTOR, "band_factor": list(BAND_FACTOR),
        "decision": {
            "width_effect": "exact two-sided McNemar p<0.05 and b>c (vs beam8)",
            "paired_ratio": "Wilson when c=0; otherwise 100000-draw four-cell bootstrap",
            "gap": "width effect AND CI excludes 1 AND CI excludes 1/17 -> GAP_CONFIRMED",
            "band": "CI inside [1/22.1, 1/13.08] REPRODUCED; wholly outside NOT_REPRODUCED; else INCONCLUSIVE",
            "k1": "same-shot McNemar vs frozen GB7 beam64_640iters mask; p<0.05 and b>c -> improvement",
            "verdict": "GAP_CONFIRMED -> FROZEN-NEGATIVE; 1/17 not excluded and improvement -> FROZEN-CERTIFIED; else FROZEN-INCONCLUSIVE",
        },
        "beam_threads": threads,
        "pre_statement": "GB9",
    }
    campaign = Campaign(config, adopt=run_dir)
    cdir = campaign.dir
    print(f"campaign: {cdir}", flush=True)
    snapshots = (
        (Path(__file__), "gateb_gb9_32res.py"),
        (Path(__file__).resolve().parent / "gateb_gb7_paired.py", "gateb_gb7_paired.py"),
        (Path(__file__).resolve().parent / "gateb_gb5a_ladder.py", "gateb_gb5a_ladder.py"),
        (Path(__file__).resolve().parent / "gateb_gb9_sizing.py", "gateb_gb9_sizing.py"),
        (BOOTSTRAP_SOURCE, "source_bootstrap.py"),
        (DEM_MATRICES_SOURCE, "source_dem_matrices.py"),
        (BEAM_SEARCH_SOURCE, "source_beam_search.py"),
        (BEAM_SOURCE, "source_beam8.cpp"),
        (TARGET / "pre_statement.md", "pre_statement.md"),
    )
    for source, name in snapshots:
        (cdir / name).write_bytes(source.read_bytes())

    circuit = load_circuit(P, BASIS)
    dem = load_dem(P, BASIS)
    if (dem.num_detectors, dem.num_errors, dem.num_observables) != (936, 8784, 12):
        raise AssertionError("wrong pinned p=1e-3 Z DEM shape")
    (cdir / "pinned.dem").write_text(str(dem))
    _h, _a, lam = dem_to_matrices(dem, merge=False)
    np.save(cdir / "pinned_lam.npy", lam)
    campaign.write_manifest({
        "circuit_sha256": circuit_sha(P, BASIS),
        "dem_sha256": dem_sha(P, BASIS),
        "dem_dimensions": {"detectors": dem.num_detectors, "errors": dem.num_errors,
                           "observables": dem.num_observables},
        "runner_sha256": sha256_file(cdir / "gateb_gb9_32res.py"),
        "sizing_sha256": sha256_file(SIZING),
        "beam_source_sha256": sha256_file(BEAM_SOURCE),
        "beam_binary": BEAM_NR_BINARY.name,
        "beam_binary_sha256": sha256_file(BEAM_NR_BINARY),
        "frozen_k1_binary_sha256": sha256_file(BEAM_DIR / "beam8_cpp"),
        "beam_equivalence_evidence": equivalence,
        "prefix_campaign": gb7,
        "versions": {"numpy": np.__version__, "stim": stim.__version__,
                     "ldpc": importlib.metadata.version("ldpc")},
        "pre_statement": "pre_statement.md Revision GB9",
    })

    # ---- shared stream: resample, must equal the frozen GB7 stream byte-for-byte
    shots_path = cdir / "shots_bposd.bin"
    inv_path = cdir / "sample_inventory.json"
    frozen_shots_sha = gb7["verified_files"]["shots_bposd.bin"]
    if inv_path.is_file() and shots_path.is_file():
        inventory = json.loads(inv_path.read_text())
        if sha256_file(shots_path) != inventory["shared_stream"]["sha256"]:
            raise ValueError("shots_bposd.bin drifted from its inventory")
        shared = inventory["shared_stream"]
        print("shared stream: reused (hash verified)", flush=True)
    else:
        print(f"sampling {SHOTS} frozen shared shots ...", flush=True)
        shared = write_shots(shots_path, circuit, SHOTS, seed)
        inventory = {"sampling_frozen_before_decode": True, "shared_stream": shared,
                     "identical_to_frozen_gb7_stream": shared["sha256"] == frozen_shots_sha}
        inv_path.write_text(json.dumps(inventory, indent=2))
    if shared["sha256"] != frozen_shots_sha:
        raise AssertionError("GB9 shared stream differs from the frozen GB7 stream: VOID")
    print("stream identity with frozen GB7: PASS", flush=True)

    shards = write_shards(cdir, "shots_bposd.bin", nshards)
    nobs = dem.num_observables

    # ---- beam8 (K=1 regression on 1e8 real shots, new binary)
    print("beam8: decoding shards ...", flush=True)
    beam8 = decode_arm(cdir, "beam8", shards, threads, nobs)
    frozen_beam8_sha = gb7["verified_files"]["beam8_predictions.bin"]
    if beam8["predictions_sha256"] != frozen_beam8_sha:
        raise AssertionError("beam8 predictions differ from frozen GB7 (new binary K=1 path): VOID")
    beam8_mask, beam8_per_obs = failure_mask_from_predictions(shots_path, cdir / "beam8_predictions.bin")
    frozen_beam8_mask = read_failmask(GB7_CAMPAIGN / "beam8_failmask.bin", SHOTS)
    if not np.array_equal(beam8_mask, frozen_beam8_mask):
        raise AssertionError("beam8 failure mask differs from frozen GB7: VOID")
    beam8_mask_sha = write_failmask(cdir / "beam8_failmask.bin", beam8_mask)
    campaign.append_result({
        "arm": "beam8-cpp-bit-exact", "rung": "beam8", "shots": SHOTS,
        "failures": int(beam8_mask.sum()),
        "raw_logical_failure_rate": float(beam8_mask.sum()) / SHOTS,
        "raw_ler_ci95_wilson": list(wilson_interval(int(beam8_mask.sum()), SHOTS)),
        "observable_mismatch_counts": beam8_per_obs.astype(int).tolist(),
        "decode": beam8, "failure_mask": "beam8_failmask.bin", "failure_mask_sha256": beam8_mask_sha,
        "identical_to_frozen_gb7_predictions": True, "identical_to_frozen_gb7_mask": True,
    })
    print(f"beam8: {int(beam8_mask.sum())} failures; identity with frozen GB7: PASS", flush=True)

    # ---- beam64_32res
    print("beam64_32res: decoding shards ...", flush=True)
    rung = decode_arm(cdir, "beam64_32res", shards, threads, nobs)
    rung_mask, rung_per_obs = failure_mask_from_predictions(shots_path, cdir / "beam64_32res_predictions.bin")
    rung_mask_sha = write_failmask(cdir / "beam64_32res_failmask.bin", rung_mask)
    campaign.append_result({
        "arm": "beam64_32res-cpp-bit-exact", "rung": "beam64_32res", "shots": SHOTS,
        "failures": int(rung_mask.sum()),
        "raw_logical_failure_rate": float(rung_mask.sum()) / SHOTS,
        "raw_ler_ci95_wilson": list(wilson_interval(int(rung_mask.sum()), SHOTS)),
        "observable_mismatch_counts": rung_per_obs.astype(int).tolist(),
        "decode": rung, "failure_mask": "beam64_32res_failmask.bin",
        "failure_mask_sha256": rung_mask_sha,
    })
    frozen_beam64_mask = read_failmask(GB7_CAMPAIGN / "beam64_failmask.bin", SHOTS)
    cells8 = paired_cells(beam8_mask, rung_mask)
    cells64 = paired_cells(frozen_beam64_mask, rung_mask)
    k1_ratio = float(gb7["beam64_decision"]["paired_ratio_rung_over_beam8"]["ratio"])
    decision = decide(cells8, cells64, SHOTS, k1_ratio)
    summary = {
        "gate": "B-gb9-beam64-32res",
        "p": P, "basis": BASIS, "shots_shared": SHOTS,
        "prefix_campaign_id": GB7_CAMPAIGN_ID,
        "stream_identity": {"shots_sha256": shared["sha256"], "identical_to_frozen_gb7": True},
        "beam8_identity": {"predictions_sha256": beam8["predictions_sha256"],
                           "identical_to_frozen_gb7": True},
        "arms": {"beam8": {"failures": int(beam8_mask.sum())},
                 "beam64_32res": {"failures": int(rung_mask.sum())}},
        "frozen_gb7_beam64_failures": int(frozen_beam64_mask.sum()),
        "decision": decision,
        "terminal_outcome": decision["outcome"],
        "band_outcome": decision["band_outcome"],
        "terminal_verdict_by_prereg_mapping": decision["terminal_verdict_by_prereg_mapping"],
    }
    campaign.close(summary)
    print(json.dumps({"campaign": str(cdir),
                      "failures": summary["arms"],
                      "outcome": decision["outcome"],
                      "band": decision["band_outcome"],
                      "vs_k1": decision["vs_beam64_640iters"]["improvement_over_k1_confirmed"],
                      "verdict": decision["terminal_verdict_by_prereg_mapping"]}, indent=2))
    return cdir


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, help="run directory minted by scripts/campaign.py init")
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--shards", type=int, default=SHOTS // SHARD)
    args = ap.parse_args()
    try:
        os.nice(2)
    except OSError:
        pass
    run(args)


if __name__ == "__main__":
    main()
