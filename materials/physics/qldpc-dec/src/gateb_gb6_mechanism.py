"""GB6 beam-expansion mechanism companion (pre_statement.md GB6 + GB5a).

Paired dense-count study at p=3e-3 basis Z on the CLEAN derived artifact
(noise-argument-only rescale; 12 observables). ONE shared seeded stream is
decoded by all four decoders, so every comparison is paired shot-by-shot:

    bp30+osd (ldpc BpOsdDecoder, pinned IonQ config, heterogeneous channel)
    beam8    (8, 30, 20, 10, num_results=1)
    beam32   (32, 40, 30, 10, num_results=1)
    beam64   (64, 40, 30, 20, num_results=1)

Frozen pre-decode clauses (GB6):
  M1 expansion monotonicity: f(beam8) >= f(beam32) >= f(beam64) AND the
     beam8/beam64 failure-count ratio CI excludes 1 from above.
  M2 sign stability: beam8-vs-BP+OSD ratio CI at 3e-3 excludes 1, and its
     direction matches the frozen p=1e-3 point ratio (0.25 < 1).
  M3 ladder sanity: f(beam64) <= f(beam32).

Where the ladder at p=1e-3 is starved of failures, this point produces dense
counts, so a real width effect must show up here if it exists at all.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gateb_gb5a_ladder import (  # noqa: E402
    BEAM_BINARY, BEAM_SOURCE, BOOTSTRAP_SOURCE, BP_OSD_SOURCE, DEM_MATRICES_SOURCE,
    BEAM_SEARCH_SOURCE, BOOTSTRAP_DRAWS, IONQ_BPOSD, RUNGS, ROUNDS,
    equivalence_evidence, failure_mask_from_predictions, mcnemar_exact,
    read_failmask, run_beam, sha256_file, shard_bounds, write_failmask,
    write_shots, load_shots_slice,
)
from qldpc_dec.bootstrap import ratio_ci, wilson_interval  # noqa: E402
from qldpc_dec.bp_osd import BpOsdBatchDecoder  # noqa: E402
from qldpc_dec.circuits import circuit_sha, dem_sha, load_circuit, load_dem  # noqa: E402
from qldpc_dec.dem_matrices import dem_to_matrices  # noqa: E402
from qldpc_dec.runner import Campaign  # noqa: E402
from qldpc_dec.seeds import derive_seed, sampling_seed  # noqa: E402

P = 3e-3
BASIS = "Z"
BASE_SEED = 20260829
SHARED_ARM = "gb6-shared"
GB5_BEAM8_POINT_RATIO_AT_1EM3 = 0.25  # frozen GB3 campaign point ratio


def decode_bposd_shard_3em3(
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
    det, obs = load_shots_slice(campaign_dir / "shots_shared.bin", start, count)
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


def paired_ratio(num_f: int, den_f: int, n: int, seed: int) -> dict[str, object]:
    """Ratio of two failure counts on the SAME n shots, with GB2-semantics CI."""
    if den_f == 0 or num_f == 0:
        return {"ratio": None if den_f == 0 else 0.0, "ci95": [0.0, None],
                "excludes_one": False, "numerator_failures": num_f,
                "denominator_failures": den_f}
    lo, hi = ratio_ci(num_f, n, den_f, n, n_boot=BOOTSTRAP_DRAWS,
                      rng=np.random.default_rng(seed))
    excludes_one = bool(lo > 1.0 or (np.isfinite(hi) and hi < 1.0))
    return {
        "ratio": float(num_f / den_f),
        "ci95": [float(lo), float(hi) if np.isfinite(hi) else None],
        "ci95_upper_unbounded": not bool(np.isfinite(hi)),
        "excludes_one": excludes_one,
        "numerator_failures": num_f, "denominator_failures": den_f,
        "bootstrap_seed": seed,
    }


def run(args: argparse.Namespace) -> Path:
    shots = int(args.shots)
    shards = int(args.shards)
    workers = int(args.workers)
    threads = int(args.threads)
    evidence = equivalence_evidence()

    config = {
        "run": "gateB-gb6-mechanism" + ("-smoke" if args.smoke else ""),
        "p": P, "basis": BASIS, "rounds": ROUNDS, "base_seed": BASE_SEED,
        "shots_shared": shots,
        "design": "one shared seeded stream decoded by all four decoders (paired)",
        "shared_arm": SHARED_ARM,
        "rungs": {r: {**{k: RUNGS[r][k] for k in ("beam_width", "initial_iters",
                                                  "iters_per_round", "max_rounds")},
                      "num_results": 1} for r in RUNGS},
        "bposd": {**dict(IONQ_BPOSD),
                  "prior_model": "merged DEM heterogeneous error_channel",
                  "shards": shards, "workers": workers},
        "clauses": {
            "M1": "f(beam8) >= f(beam32) >= f(beam64) and beam8/beam64 ratio CI excludes 1",
            "M2": "beam8/bp30+osd ratio CI excludes 1 with the same direction as the frozen 1e-3 point (<1)",
            "M3": "f(beam64) <= f(beam32)",
        },
        "bootstrap": {"draws": BOOTSTRAP_DRAWS},
        "pre_statement": "GB6 + GB5a",
    }
    campaign = Campaign(config, adopt=Path(args.run_dir))
    cdir = campaign.dir
    print(f"campaign: {cdir}", flush=True)
    for src, dst in ((Path(__file__), "gateb_gb6_mechanism.py"),
                     (Path(__file__).resolve().parent / "gateb_gb5a_ladder.py",
                      "gateb_gb5a_ladder.py"),
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
        raise AssertionError(
            f"wrong 3e-3 DEM shape {(dem.num_detectors, dem.num_errors, dem.num_observables)}; "
            "the clean argrescale artifact must be used")
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
    }

    seed = sampling_seed(BASE_SEED, P, BASIS, SHARED_ARM)

    campaign.write_manifest({
        "circuit_sha256": circuit_sha(P, BASIS), "dem_sha256": dem_sha(P, BASIS),
        "circuit_file": "BB_144_144_12_memory_Z_p0.003_sr12_derived_p1e-3_argrescale.stim",
        "circuit_derivation": "noise-argument-only rescale (0.001)->(0.003) of the pinned IonQ Z p=1e-3 circuit",
        "dem_dimensions": {"detectors": dem.num_detectors, "errors": dem.num_errors,
                           "observables": dem.num_observables},
        "sampling_seed": seed,
        "runner_sha256": sha256_file(cdir / "gateb_gb6_mechanism.py"),
        "ladder_runner_sha256": sha256_file(cdir / "gateb_gb5a_ladder.py"),
        "bootstrap_source_sha256": sha256_file(cdir / "source_bootstrap.py"),
        "bp_osd_source_sha256": sha256_file(cdir / "source_bp_osd.py"),
        "beam_source_sha256": sha256_file(BEAM_SOURCE),
        "beam_binary_sha256": sha256_file(BEAM_BINARY),
        "beam_equivalence_evidence": evidence,
        "bposd_prior_check": prior_check,
        "versions": {"numpy": np.__version__, "stim": stim.__version__,
                     "ldpc": importlib.metadata.version("ldpc")},
        "pre_statement": "pre_statement.md Gate B revisions GB6 + GB5a",
    })
    shared = write_shots(cdir / "shots_shared.bin", circuit, shots, seed)
    (cdir / "sample_inventory.json").write_text(json.dumps(
        {"sampling_frozen_before_decode": True, "shared_stream": shared}, indent=2))

    # ---- BP+OSD on the shared stream ----
    print(f"BP+OSD on shared stream: {shots} shots / {shards} shards ...", flush=True)
    bounds = shard_bounds(shots, shards)
    ctx = multiprocessing.get_context("spawn")
    results: list[dict[str, object]] = []
    t0 = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as pool:
        futs = [pool.submit(decode_bposd_shard_3em3, str(cdir), s, st, c, str(shared["sha256"]))
                for s, st, c in bounds]
        for fut in concurrent.futures.as_completed(futs):
            results.append(fut.result())
    bp_wall = time.perf_counter() - t0
    results.sort(key=lambda r: int(r["shard"]))
    bp_failures = sum(int(r["failures"]) for r in results)
    bp_mask = np.concatenate([read_failmask(cdir / str(r["failure_mask"]), int(r["shots"]))
                              for r in results])
    if bp_mask.size != shots or int(bp_mask.sum()) != bp_failures:
        raise AssertionError("stitched BP+OSD mask disagrees with shard counts")
    write_failmask(cdir / "bposd_failmask.bin", bp_mask)
    core_s = sum(float(r["wall_s"]) for r in results)
    bposd = {
        "arm": "bp30+osd", "shots": shots, "failures": bp_failures,
        "raw_logical_failure_rate": bp_failures / shots,
        "raw_ler_ci95_wilson": list(wilson_interval(bp_failures, shots)),
        "per_round_ler": bp_failures / shots / ROUNDS,
        "decode_core_s": core_s, "decode_ms_per_shot": 1000.0 * core_s / shots,
        "wall_s": bp_wall, "shards": len(results), "decoder": dict(IONQ_BPOSD),
    }
    campaign.append_result(bposd)
    print(f"  BP+OSD: {bp_failures}/{shots}", flush=True)

    # ---- beam rungs on the same shared stream ----
    masks: dict[str, np.ndarray] = {"bp30+osd": bp_mask}
    arms: dict[str, dict[str, object]] = {"bp30+osd": bposd}
    for rung in RUNGS:
        print(f"{rung} on shared stream ...", flush=True)
        meta = run_beam(
            cdir, rung, "shots_shared.bin", f"{rung}_predictions.bin", rung,
            threads, str(shared["sha256"]))
        failed, per_obs = failure_mask_from_predictions(
            cdir / "shots_shared.bin", cdir / f"{rung}_predictions.bin")
        mask_sha = write_failmask(cdir / f"{rung}_failmask.bin", failed)
        f = int(failed.sum())
        arm = {
            "arm": f"{rung}-cpp-bit-exact", "rung": rung, "shots": shots, "failures": f,
            "raw_logical_failure_rate": f / shots,
            "raw_ler_ci95_wilson": list(wilson_interval(f, shots)),
            "per_round_ler": f / shots / ROUNDS,
            "observable_mismatch_counts": per_obs.astype(int).tolist(),
            "decoder_params": {k: RUNGS[rung][k] for k in
                               ("beam_width", "initial_iters", "iters_per_round", "max_rounds")},
            "decode": meta, "failure_mask": f"{rung}_failmask.bin",
            "failure_mask_sha256": mask_sha,
        }
        campaign.append_result(arm)
        masks[rung] = failed
        arms[rung] = arm
        print(f"  {rung}: {f}/{shots}", flush=True)

    names = list(arms)
    discordance = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            fa, fb = masks[a], masks[b]
            n_ab, n_ba = int((fa & ~fb).sum()), int((fb & ~fa).sum())
            discordance[f"{a}|{b}"] = {
                "a_fails_b_ok": n_ab, "b_fails_a_ok": n_ba,
                "both_fail": int((fa & fb).sum()),
                "mcnemar_exact_p": mcnemar_exact(n_ab, n_ba),
            }
    identity = {}
    rung_list = list(RUNGS)
    for i, a in enumerate(rung_list):
        for b in rung_list[i + 1:]:
            sa = sha256_file(cdir / f"{a}_predictions.bin")
            sb = sha256_file(cdir / f"{b}_predictions.bin")
            identity[f"{a}|{b}"] = {"identical_predictions": sa == sb, "sha256_a": sa, "sha256_b": sb}

    f8 = int(arms["beam8"]["failures"])
    f32 = int(arms["beam32"]["failures"])
    f64 = int(arms["beam64"]["failures"])
    r_8_over_64 = paired_ratio(f8, f64, shots, derive_seed(BASE_SEED, "gb6-ratio", "beam8/beam64"))
    r_8_over_bp = paired_ratio(f8, bp_failures, shots, derive_seed(BASE_SEED, "gb6-ratio", "beam8/bposd"))
    r_32_over_bp = paired_ratio(f32, bp_failures, shots, derive_seed(BASE_SEED, "gb6-ratio", "beam32/bposd"))
    r_64_over_bp = paired_ratio(f64, bp_failures, shots, derive_seed(BASE_SEED, "gb6-ratio", "beam64/bposd"))

    m1 = bool(f8 >= f32 >= f64) and bool(r_8_over_64["excludes_one"]) and bool(
        (r_8_over_64["ratio"] or 0) > 1.0)
    direction_matches = bool(
        r_8_over_bp["ratio"] is not None
        and ((r_8_over_bp["ratio"] < 1.0) == (GB5_BEAM8_POINT_RATIO_AT_1EM3 < 1.0)))
    m2 = bool(r_8_over_bp["excludes_one"]) and direction_matches
    m3 = bool(f64 <= f32)
    clauses = {
        "M1_expansion_monotonicity": {"pass": m1, "ordering_ok": bool(f8 >= f32 >= f64),
                                      "beam8_over_beam64": r_8_over_64},
        "M2_sign_stability": {"pass": m2, "beam8_over_bposd": r_8_over_bp,
                              "frozen_1e-3_point_ratio": GB5_BEAM8_POINT_RATIO_AT_1EM3,
                              "direction_matches": direction_matches},
        "M3_ladder_sanity": {"pass": m3, "beam64_failures": f64, "beam32_failures": f32},
    }
    summary = {
        "gate": "B-gb6-mechanism-paired",
        "p": P, "basis": BASIS, "shots_shared": shots,
        "arms": arms,
        "failures": {k: int(v.sum()) for k, v in masks.items()},
        "paired_discordance": discordance,
        "prediction_identity": identity,
        "ratios_vs_bposd": {"beam8": r_8_over_bp, "beam32": r_32_over_bp, "beam64": r_64_over_bp},
        "clauses": clauses,
        "mechanism_verdict": (
            "EXPANSION_EFFECT_PRESENT" if (m1 and m3) else
            "NO_MEASURABLE_EXPANSION_EFFECT" if all(
                discordance[k]["a_fails_b_ok"] == 0 and discordance[k]["b_fails_a_ok"] == 0
                for k in ("beam8|beam32", "beam8|beam64", "beam32|beam64")) else
            "PARTIAL_OR_INCONSISTENT_EXPANSION_EFFECT"),
        "scope": "mechanism companion; no GB5 band decision is altered by these clauses",
    }
    campaign.close(summary)
    print(json.dumps({"campaign": str(cdir), "failures": summary["failures"],
                      "verdict": summary["mechanism_verdict"],
                      "clauses": {k: v["pass"] for k, v in clauses.items()}}, indent=2))
    return cdir


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shots", type=int, default=200_000)
    ap.add_argument("--shards", type=int, default=48)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--threads", type=int, default=26)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--run-dir", required=True,
                    help="run directory minted by scripts/campaign.py init")
    args = ap.parse_args()
    if args.smoke:
        if args.shots == 200_000:
            args.shots = 2_000
        if args.shards == 48:
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
