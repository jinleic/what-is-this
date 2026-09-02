"""GB7 paired beam-expansion escalation (pre_statement.md Revision GB7).

Decode beam8, beam32, and beam64 on the first 1e8 shots of the frozen
bp30+osd stream at p=1e-3, basis Z. The first 2e7 shots and beam predictions
must match the frozen GB5a paired campaign byte-for-byte. Width significance
is decided by exact McNemar; paired-ratio intervals determine whether the
measured factor is quantitatively incompatible with the published factor.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import sys
from pathlib import Path
from types import MappingProxyType

import numpy as np
import stim

TARGET = Path(__file__).resolve().parents[1]
SRC = TARGET / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gateb_gb5a_ladder import (  # noqa: E402
    BASE_SEED,
    BASIS,
    BEAM_BINARY,
    BEAM_SEARCH_SOURCE,
    BEAM_SOURCE,
    BOOTSTRAP_SOURCE,
    CHUNK,
    DEM_MATRICES_SOURCE,
    RUNGS,
    check_predictions_header,
    equivalence_evidence,
    failure_mask_from_predictions,
    mcnemar_exact,
    read_shots_header,
    run_beam,
    sha256_file,
    write_failmask,
    write_shots,
)
from qldpc_dec.bootstrap import wilson_interval  # noqa: E402
from qldpc_dec.circuits import circuit_sha, dem_sha, load_circuit, load_dem  # noqa: E402
from qldpc_dec.dem_matrices import dem_to_matrices  # noqa: E402
from qldpc_dec.runner import Campaign  # noqa: E402
from qldpc_dec.seeds import derive_seed, sampling_seed  # noqa: E402

P = 1e-3
ROUNDS = 12
SHOTS = 100_000_000
PREFIX_SHOTS = 20_000_000
BOOTSTRAP_DRAWS = 100_000
MCNEMAR_ALPHA = 0.05
PREFIX_CAMPAIGN_ID = "20260901T040018Z_55e4ac01_f81265905380"
PREFIX_CAMPAIGN = TARGET / "campaigns" / PREFIX_CAMPAIGN_ID
PUBLISHED_RECIPROCALS = MappingProxyType({
    "beam32": 1.0 / 5.6,
    # Revision GB8: mis-attributed target. The configuration run
    # (beam64_640iters, num_results=1) is published at 7.0x; 17x is
    # beam64_32res_640iters (num_results=32). Kept verbatim because the frozen
    # GB7 decision and verify_gateb_gb5.py's replay depend on it; the
    # corrected re-target lives in gateb_gb8_retarget.py.
    "beam64": 1.0 / 17.0,
})


def sha256_region(path: Path, offset: int, size: int) -> str:
    """Hash exactly one byte region, rejecting truncation."""
    if offset < 0 or size < 0:
        raise ValueError("region offset and size must be non-negative")
    digest = hashlib.sha256()
    remaining = size
    with path.open("rb") as handle:
        handle.seek(offset)
        while remaining:
            block = handle.read(min(1 << 22, remaining))
            if not block:
                raise ValueError(f"truncated region in {path}: {remaining} bytes missing")
            digest.update(block)
            remaining -= len(block)
    return digest.hexdigest()


def frozen_checksums(run_dir: Path) -> dict[str, str]:
    path = run_dir / "sha256s.txt"
    if not path.is_file():
        raise FileNotFoundError(f"prefix campaign is not frozen: {path} missing")
    checksums: dict[str, str] = {}
    for line in path.read_text().splitlines():
        digest, separator, relative = line.partition("  ")
        if not separator or len(digest) != 64 or relative in checksums:
            raise ValueError(f"malformed frozen checksum line: {line!r}")
        checksums[relative] = digest
    return checksums


def verify_frozen_file(run_dir: Path, relative: str, checksums: dict[str, str]) -> str:
    if relative not in checksums:
        raise ValueError(f"{relative} is absent from frozen prefix checksums")
    path = run_dir / relative
    actual = sha256_file(path)
    if actual != checksums[relative]:
        raise ValueError(f"frozen prefix artifact changed: {path}")
    return actual


def validate_prefix_campaign(run_dir: Path) -> dict[str, object]:
    if run_dir.resolve() != PREFIX_CAMPAIGN.resolve():
        raise ValueError(
            f"GB7 is pinned to prefix campaign {PREFIX_CAMPAIGN}, got {run_dir}")
    status_path = run_dir / "status.json"
    if not status_path.is_file():
        raise FileNotFoundError(f"prefix campaign has no terminal status: {status_path}")
    status = json.loads(status_path.read_text())
    verdict = status.get("verdict")
    if not isinstance(verdict, str) or not verdict.startswith("FROZEN-"):
        raise ValueError(f"prefix campaign is not evidentially frozen: {status}")

    checksums = frozen_checksums(run_dir)
    required = [
        "manifest.json",
        "summary.json",
        "shots_bposd.bin",
        "paired_beam8_predictions.bin",
        "paired_beam32_predictions.bin",
        "paired_beam64_predictions.bin",
        "paired_beam8_failmask.bin",
        "paired_beam32_failmask.bin",
        "paired_beam64_failmask.bin",
    ]
    verified = {name: verify_frozen_file(run_dir, name, checksums) for name in required}
    summary = json.loads((run_dir / "summary.json").read_text())
    if int(summary.get("shots_per_arm", -1)) != PREFIX_SHOTS:
        raise ValueError("GB5a prefix campaign does not contain exactly 2e7 shots")
    paired = summary.get("paired_instrument", {})
    if int(paired.get("shots", -1)) != PREFIX_SHOTS:
        raise ValueError("GB5a paired instrument does not contain exactly 2e7 shots")
    return {
        "run_id": run_dir.name,
        "terminal_verdict": verdict,
        "sha256s_txt_sha256": sha256_file(run_dir / "sha256s.txt"),
        "verified_files": verified,
        "paired_failures": paired.get("failures"),
    }


def compare_shot_prefix(current: Path, frozen: Path) -> dict[str, object]:
    c_ndet, c_nobs, c_shots, c_det_bytes, c_obs_bytes = read_shots_header(current)
    f_ndet, f_nobs, f_shots, f_det_bytes, f_obs_bytes = read_shots_header(frozen)
    if (c_ndet, c_nobs) != (f_ndet, f_nobs):
        raise AssertionError("GB7 and GB5a shot dimensions differ")
    if c_shots != SHOTS or f_shots != PREFIX_SHOTS:
        raise AssertionError(f"wrong shot counts for prefix check: {c_shots}, {f_shots}")
    if (c_det_bytes, c_obs_bytes) != (f_det_bytes, f_obs_bytes):
        raise AssertionError("GB7 and GB5a packed row widths differ")

    regions = {
        "detectors": (
            sha256_region(current, 24, PREFIX_SHOTS * c_det_bytes),
            sha256_region(frozen, 24, PREFIX_SHOTS * f_det_bytes),
        ),
        "observables": (
            sha256_region(
                current,
                24 + c_shots * c_det_bytes,
                PREFIX_SHOTS * c_obs_bytes,
            ),
            sha256_region(
                frozen,
                24 + f_shots * f_det_bytes,
                PREFIX_SHOTS * f_obs_bytes,
            ),
        ),
    }
    result: dict[str, object] = {}
    for name, (current_sha, frozen_sha) in regions.items():
        if current_sha != frozen_sha:
            raise AssertionError(f"GB7 {name} prefix differs from frozen GB5a")
        result[name] = {
            "bytes": PREFIX_SHOTS * (c_det_bytes if name == "detectors" else c_obs_bytes),
            "sha256": current_sha,
            "byte_match": True,
        }
    return result


def compare_prediction_prefix(current: Path, frozen: Path, nobs: int) -> dict[str, object]:
    check_predictions_header(current, nobs, SHOTS)
    check_predictions_header(frozen, nobs, PREFIX_SHOTS)
    row_bytes = ((nobs + 63) // 64) * 8
    size = PREFIX_SHOTS * row_bytes
    current_sha = sha256_region(current, 12, size)
    frozen_sha = sha256_region(frozen, 12, size)
    if current_sha != frozen_sha:
        raise AssertionError(f"GB7 prediction prefix differs from frozen GB5a: {current.name}")
    return {"bytes": size, "sha256": current_sha, "byte_match": True}


def paired_cells(reference: np.ndarray, candidate: np.ndarray) -> dict[str, int]:
    if reference.dtype != np.bool_ or candidate.dtype != np.bool_:
        raise TypeError("paired failure masks must be boolean")
    if reference.shape != candidate.shape:
        raise ValueError("paired failure masks must have identical shape")
    both = int((reference & candidate).sum())
    ref_only = int((reference & ~candidate).sum())
    candidate_only = int((~reference & candidate).sum())
    neither = int(reference.size - both - ref_only - candidate_only)
    return {
        "both_fail": both,
        "reference_fails_candidate_ok": ref_only,
        "reference_ok_candidate_fails": candidate_only,
        "neither_fails": neither,
    }


def paired_ratio_interval(cells: dict[str, int], rung: str, shots: int) -> dict[str, object]:
    both = int(cells["both_fail"])
    ref_only = int(cells["reference_fails_candidate_ok"])
    candidate_only = int(cells["reference_ok_candidate_fails"])
    neither = int(cells["neither_fails"])
    counts = np.asarray([both, ref_only, candidate_only, neither], dtype=np.int64)
    if np.any(counts < 0) or int(counts.sum()) != shots:
        raise ValueError(f"invalid paired table for {rung}: {counts.tolist()}")
    f8 = both + ref_only
    fr = both + candidate_only
    seed = derive_seed(BASE_SEED, "gb7-paired-ratio", rung, shots)
    if f8 == 0:
        return {
            "method": "undefined_no_beam8_failures",
            "ratio": None,
            "ci95": [None, None],
            "beam8_failures": 0,
            "rung_failures": fr,
            "bootstrap_seed": seed,
        }

    ratio = fr / f8
    if candidate_only == 0:
        lo, hi = wilson_interval(both, f8)
        method = "Wilson on both_fail / beam8_failures (exact subset case c=0)"
        draws = None
        valid_draws = None
    else:
        rng = np.random.default_rng(seed)
        samples = rng.multinomial(shots, counts / shots, size=BOOTSTRAP_DRAWS)
        denominators = samples[:, 0] + samples[:, 1]
        numerators = samples[:, 0] + samples[:, 2]
        valid = denominators > 0
        boot = numerators[valid] / denominators[valid]
        if boot.size == 0:
            raise RuntimeError(f"all paired bootstrap draws had zero beam8 failures for {rung}")
        lo, hi = np.quantile(boot, [0.025, 0.975])
        method = "four-cell paired multinomial percentile bootstrap"
        draws = BOOTSTRAP_DRAWS
        valid_draws = int(boot.size)
    return {
        "method": method,
        "ratio": float(ratio),
        "ci95": [float(lo), float(hi)],
        "beam8_failures": f8,
        "rung_failures": fr,
        "bootstrap_draws": draws,
        "valid_bootstrap_draws": valid_draws,
        "bootstrap_seed": seed,
    }


def rung_decision(rung: str, cells: dict[str, int], shots: int) -> dict[str, object]:
    interval = paired_ratio_interval(cells, rung, shots)
    b = int(cells["reference_fails_candidate_ok"])
    c = int(cells["reference_ok_candidate_fails"])
    p_value = mcnemar_exact(b, c)
    width_effect = bool(p_value < MCNEMAR_ALPHA and b > c)
    ratio = interval["ratio"]
    lo, hi = interval["ci95"]
    reciprocal = float(PUBLISHED_RECIPROCALS[rung])
    if ratio is None or lo is None or hi is None:
        excludes_one = False
        published_excluded = False
        factor = None
        factor_ci: list[float | None] = [None, None]
    else:
        excludes_one = bool(hi < 1.0 or lo > 1.0)
        published_excluded = bool(reciprocal < lo or reciprocal > hi)
        factor = None if ratio == 0 else float(1.0 / ratio)
        factor_ci = [float(1.0 / hi), None if lo == 0 else float(1.0 / lo)]
    gap = bool(width_effect and excludes_one and published_excluded)
    satisfied = [
        label for label, condition in (
            ("WIDTH_EFFECT_CONFIRMED", width_effect),
            ("INTERVAL_EXCLUDES_ONE", excludes_one),
            ("PUBLISHED_FACTOR_EXCLUDED", published_excluded),
        )
        if condition
    ]
    outcome = (
        "GAP_CONFIRMED" if gap else
        "NO_GAP_" + "_AND_".join(satisfied) if satisfied else
        "UNRESOLVED_AT_N"
    )
    return {
        "rung": rung,
        "paired_cells": cells,
        "mcnemar": {
            "exact_two_sided_p": p_value,
            "alpha": MCNEMAR_ALPHA,
            "directional_benefit": b > c,
            "width_effect_confirmed": width_effect,
            "b_beam8_fails_rung_ok": b,
            "c_beam8_ok_rung_fails": c,
        },
        "paired_ratio_rung_over_beam8": interval,
        "interval_excludes_one": excludes_one,
        "published_reciprocal": reciprocal,
        "published_factor": float(1.0 / reciprocal),
        "published_factor_excluded": published_excluded,
        "measured_factor_beam8_over_rung": factor,
        "measured_factor_ci95": factor_ci,
        "gap_confirmed": gap,
        "outcome": outcome,
    }


def run(args: argparse.Namespace) -> Path:
    run_dir = Path(args.run_dir).resolve()
    prefix_dir = Path(args.prefix_run).resolve()
    threads = int(args.threads)
    if threads <= 0:
        raise ValueError("threads must be positive")
    # write_shots draws in CHUNK-sized sampler calls, so the GB7 stream repeats
    # the GB5a call sequence verbatim only when both counts are chunk-aligned.
    # Unaligned counts silently draw a different stim stream; assert instead.
    if SHOTS % CHUNK or PREFIX_SHOTS % CHUNK:
        raise ValueError(
            f"shot counts must be multiples of CHUNK={CHUNK} for prefix "
            f"continuity; got SHOTS={SHOTS}, PREFIX_SHOTS={PREFIX_SHOTS}")

    prefix_evidence = validate_prefix_campaign(prefix_dir)
    equivalence = equivalence_evidence()
    seed = sampling_seed(BASE_SEED, P, BASIS, "bp30+osd")
    config = {
        "run": "gateB-gb7-paired-expansion",
        "p": P,
        "basis": BASIS,
        "rounds": ROUNDS,
        "base_seed": BASE_SEED,
        "sampling_seed": seed,
        "shots_shared": SHOTS,
        "prefix_shots": PREFIX_SHOTS,
        "prefix_campaign_id": PREFIX_CAMPAIGN_ID,
        "rungs": {
            rung: {
                **{key: RUNGS[rung][key] for key in (
                    "beam_width", "initial_iters", "iters_per_round", "max_rounds")},
                "num_results": 1,
            }
            for rung in RUNGS
        },
        "decision": {
            "width_effect": "exact two-sided McNemar p<0.05 and b>c",
            "paired_ratio": "Wilson when c=0; otherwise 100000-draw four-cell bootstrap",
            "gap": "McNemar width effect and paired CI excludes both 1 and published reciprocal",
        },
        "beam_threads": threads,
        "pre_statement": "GB7",
    }
    campaign = Campaign(config, adopt=run_dir)
    cdir = campaign.dir
    print(f"campaign: {cdir}", flush=True)

    snapshots = (
        (Path(__file__), "gateb_gb7_paired.py"),
        (Path(__file__).resolve().parent / "gateb_gb5a_ladder.py", "gateb_gb5a_ladder.py"),
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
        "dem_dimensions": {
            "detectors": dem.num_detectors,
            "errors": dem.num_errors,
            "observables": dem.num_observables,
        },
        "runner_sha256": sha256_file(cdir / "gateb_gb7_paired.py"),
        "ladder_runner_sha256": sha256_file(cdir / "gateb_gb5a_ladder.py"),
        "beam_source_sha256": sha256_file(BEAM_SOURCE),
        "beam_binary_sha256": sha256_file(BEAM_BINARY),
        "beam_equivalence_evidence": equivalence,
        "prefix_campaign": prefix_evidence,
        "versions": {
            "numpy": np.__version__,
            "stim": stim.__version__,
            "ldpc": importlib.metadata.version("ldpc"),
        },
        "pre_statement": "pre_statement.md Revision GB7",
    })

    print(f"sampling {SHOTS} frozen shared shots ...", flush=True)
    shared = write_shots(cdir / "shots_bposd.bin", circuit, SHOTS, seed)
    shot_prefix = compare_shot_prefix(
        cdir / "shots_bposd.bin", prefix_dir / "shots_bposd.bin")
    sample_inventory = {
        "sampling_frozen_before_decode": True,
        "shared_stream": shared,
        "prefix_continuity": shot_prefix,
    }
    (cdir / "sample_inventory.json").write_text(json.dumps(sample_inventory, indent=2))
    print("shot prefix continuity: PASS", flush=True)

    masks: dict[str, np.ndarray] = {}
    arms: dict[str, dict[str, object]] = {}
    prediction_prefix: dict[str, object] = {}
    for rung in RUNGS:
        print(f"{rung}: decoding {SHOTS} shared shots ...", flush=True)
        prediction_name = f"{rung}_predictions.bin"
        meta = run_beam(
            cdir, rung, "shots_bposd.bin", prediction_name, rung, threads,
            str(shared["sha256"]))
        prefix = compare_prediction_prefix(
            cdir / prediction_name,
            prefix_dir / f"paired_{rung}_predictions.bin",
            dem.num_observables,
        )
        failed, per_observable = failure_mask_from_predictions(
            cdir / "shots_bposd.bin", cdir / prediction_name)
        frozen_prefix_mask = np.frombuffer(
            (prefix_dir / f"paired_{rung}_failmask.bin").read_bytes(), dtype=np.uint8)
        frozen_prefix_mask = np.unpackbits(
            frozen_prefix_mask, bitorder="little")[:PREFIX_SHOTS].astype(bool)
        if not np.array_equal(failed[:PREFIX_SHOTS], frozen_prefix_mask):
            raise AssertionError(f"GB7 {rung} failure-mask prefix differs from frozen GB5a")
        mask_name = f"{rung}_failmask.bin"
        mask_sha = write_failmask(cdir / mask_name, failed)
        failures = int(failed.sum())
        arm = {
            "arm": f"{rung}-cpp-bit-exact",
            "rung": rung,
            "shots": SHOTS,
            "failures": failures,
            "raw_logical_failure_rate": failures / SHOTS,
            "raw_ler_ci95_wilson": list(wilson_interval(failures, SHOTS)),
            "per_round_ler": failures / SHOTS / ROUNDS,
            "observable_mismatch_counts": per_observable.astype(int).tolist(),
            "decode": meta,
            "failure_mask": mask_name,
            "failure_mask_sha256": mask_sha,
            "prefix_failures": int(failed[:PREFIX_SHOTS].sum()),
            "prefix_failure_mask_match": True,
        }
        campaign.append_result(arm)
        masks[rung] = failed
        arms[rung] = arm
        prediction_prefix[rung] = prefix
        print(
            f"  {rung}: {failures}/{SHOTS}; prediction and failure-mask prefix PASS",
            flush=True,
        )

    discordance: dict[str, object] = {}
    names = list(RUNGS)
    for index, first in enumerate(names):
        for second in names[index + 1:]:
            cells = paired_cells(masks[first], masks[second])
            b = int(cells["reference_fails_candidate_ok"])
            c = int(cells["reference_ok_candidate_fails"])
            discordance[f"{first}|{second}"] = {
                "reference": first,
                "candidate": second,
                **cells,
                "mcnemar_exact_two_sided_p": mcnemar_exact(b, c),
            }

    decisions = {
        rung: rung_decision(rung, paired_cells(masks["beam8"], masks[rung]), SHOTS)
        for rung in PUBLISHED_RECIPROCALS
    }
    summary = {
        "gate": "B-gb7-paired-expansion",
        "p": P,
        "basis": BASIS,
        "shots_shared": SHOTS,
        "scope": "paired beam expansion only; GB5a independent-arm band verdicts unchanged",
        "prefix_campaign_id": PREFIX_CAMPAIGN_ID,
        "prefix_continuity": {
            "shots": shot_prefix,
            "predictions": prediction_prefix,
            "failure_masks_match": {rung: True for rung in RUNGS},
        },
        "arms": arms,
        "paired_discordance": discordance,
        "decisions": decisions,
        "beam32_terminal_outcome": decisions["beam32"]["outcome"],
        "beam64_terminal_outcome": decisions["beam64"]["outcome"],
    }
    campaign.close(summary)
    print(json.dumps({
        "campaign": str(cdir),
        "failures": {rung: arm["failures"] for rung, arm in arms.items()},
        "outcomes": {rung: decision["outcome"] for rung, decision in decisions.items()},
        "prefix_continuity": "PASS",
    }, indent=2))
    return cdir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-dir",
        required=True,
        help="run directory minted by scripts/campaign.py init",
    )
    parser.add_argument(
        "--prefix-run",
        default=str(PREFIX_CAMPAIGN),
        help="frozen GB5a campaign; GB7 rejects any value except the pinned run",
    )
    parser.add_argument("--threads", type=int, default=26)
    args = parser.parse_args()
    try:
        os.nice(5)
    except OSError:
        pass
    run(args)


if __name__ == "__main__":
    main()
