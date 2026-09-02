"""Apply pre-data Gate-B Revision GB2 to a closed, GB3-valid campaign.

The source campaign remains immutable. This amendment reuses only its frozen
arm counts and original bootstrap design, correcting
positive-numerator/zero-denominator draws to retain +infinity. Campaigns that
lack the heterogeneous BP channel required by Revision GB3 are rejected:
their baseline counts are invalid, not amendable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np

from qldpc_dec.bootstrap import ratio_ci
from qldpc_dec.runner import CAMPAIGNS_ROOT, Campaign


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise TypeError(f"expected a JSON object in {path}")
    return value


def corrected_verdict(
    beam: dict[str, object],
    bposd: dict[str, object],
    *,
    draws: int,
    seed: int,
    band: tuple[float, float],
) -> dict[str, object]:
    """Recompute only the ratio interval and frozen three-way verdict."""
    beam_failures = int(beam["failures"])
    beam_shots = int(beam["shots"])
    bposd_failures = int(bposd["failures"])
    bposd_shots = int(bposd["shots"])
    if beam_shots <= 0 or bposd_shots <= 0:
        raise ValueError("both arms must contain at least one shot")

    if beam_failures == 0 or bposd_failures == 0:
        return {
            "outcome": "UNDERDETERMINED_zero_failures_in_at_least_one_arm",
            "beam8_over_bposd_ratio": None if bposd_failures == 0 else 0.0,
            "ratio_ci95": [0.0, None],
            "ratio_ci95_upper_unbounded": True,
            "target_band": list(band),
            "point_in_band": None,
            "ci_inside_band": False,
            "ci_overlaps_band": True,
            "bootstrap": {
                "method": "GB2 independent-binomial percentile",
                "draws": draws,
                "seed": seed,
                "zero_denominator_rule": (
                    "positive/zero retained as +inf; zero/zero discarded"
                ),
            },
        }

    ratio = (beam_failures * bposd_shots) / (bposd_failures * beam_shots)
    lo, hi = ratio_ci(
        beam_failures,
        beam_shots,
        bposd_failures,
        bposd_shots,
        n_boot=draws,
        rng=np.random.default_rng(seed),
    )
    point_in_band = bool(band[0] <= ratio <= band[1])
    ci_inside_band = bool(band[0] <= lo and hi <= band[1])
    ci_overlaps_band = bool(lo <= band[1] and hi >= band[0])
    if ci_inside_band:
        outcome = "REPRODUCED_ratio_CI_inside_band"
    elif not ci_overlaps_band:
        outcome = "NOT_REPRODUCED_ratio_CI_outside_band"
    else:
        outcome = "INCONCLUSIVE_ratio_CI_crosses_band"

    return {
        "outcome": outcome,
        "beam8_over_bposd_ratio": ratio,
        "ratio_ci95": [lo, float(hi) if math.isfinite(hi) else None],
        "ratio_ci95_upper_unbounded": not math.isfinite(hi),
        "target_band": list(band),
        "point_in_band": point_in_band,
        "ci_inside_band": ci_inside_band,
        "ci_overlaps_band": ci_overlaps_band,
        "bootstrap": {
            "method": "GB2 independent-binomial percentile",
            "draws": draws,
            "seed": seed,
            "zero_denominator_rule": (
                "positive/zero retained as +inf; zero/zero discarded"
            ),
        },
    }


def amend(source_campaign: Path, root: Path | None = None) -> Path:
    source_campaign = source_campaign.resolve()
    source_files = {
        name: source_campaign / name
        for name in ("manifest.json", "summary.json", "inventory.json", "results.json.gz")
    }
    missing = [str(path) for path in source_files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"source campaign is not closed: missing {missing}")

    manifest = _load_json(source_files["manifest.json"])
    summary = _load_json(source_files["summary.json"])
    inventory = _load_json(source_files["inventory.json"])
    if inventory.get("policy") != "immutable-after-close":
        raise ValueError("source campaign lacks immutable-after-close inventory")

    source_config = manifest.get("config")
    if not isinstance(source_config, dict):
        raise TypeError("source manifest config is not an object")
    expected = {
        "p": 1e-3,
        "basis": "Z",
        "rounds": 12,
        "base_seed": 20260829,
    }
    for key, value in expected.items():
        if source_config.get(key) != value:
            raise ValueError(f"source {key}={source_config.get(key)!r}, expected {value!r}")
    if source_config.get("run") not in {"gateB-pinned-1e6", "gateB-pinned-smoke"}:
        raise ValueError(f"unsupported source run {source_config.get('run')!r}")
    bposd_config = source_config.get("bposd")
    if (
        not isinstance(bposd_config, dict)
        or bposd_config.get("prior_model")
        != "merged DEM heterogeneous error_channel"
    ):
        raise ValueError(
            "source campaign is invalid under Revision GB3: "
            "missing heterogeneous BP error_channel"
        )

    beam = summary.get("beam8")
    bposd = summary.get("bposd")
    if not isinstance(beam, dict) or not isinstance(bposd, dict):
        raise TypeError("source summary lacks arm result objects")
    shots = int(source_config["shots_per_arm"])
    if int(beam.get("shots", -1)) != shots or int(bposd.get("shots", -1)) != shots:
        raise ValueError("source arm shot counts do not match the frozen manifest")

    draws = int(source_config["bootstrap_draws"])
    seed = int(source_config["bootstrap_seed"])
    band_values = source_config["target_band_beam8_over_bposd"]
    if not isinstance(band_values, list) or len(band_values) != 2:
        raise TypeError("source target band must be a two-element list")
    band = (float(band_values[0]), float(band_values[1]))

    source_sha256 = {
        name: sha256_file(path) for name, path in source_files.items()
    }
    amendment_config = {
        "run": "gateB-pinned-GB2-amendment",
        "revision": "pre_statement.md Revision GB2",
        "source_campaign": source_campaign.name,
        "source_config_hash12": manifest.get("config_hash12"),
        "source_summary_sha256": source_sha256["summary.json"],
        "shots_per_arm": shots,
        "bootstrap_draws": draws,
        "bootstrap_seed": seed,
        "target_band_beam8_over_bposd": list(band),
        "decoder_rerun": False,
        "frozen_raw_counts_reused": True,
        "source_bposd_prior_model": bposd_config["prior_model"],
    }
    campaign = Campaign(amendment_config, root=root or CAMPAIGNS_ROOT())
    implementation = Path(__file__).resolve()
    bootstrap_source = implementation.parent / "qldpc_dec" / "bootstrap.py"
    campaign.write_manifest(
        {
            "source_campaign_path": str(source_campaign),
            "source_campaign_closed_utc": inventory.get("closed_utc"),
            "source_file_sha256": source_sha256,
            "amendment_source_sha256": sha256_file(implementation),
            "bootstrap_source_sha256": sha256_file(bootstrap_source),
            "scope": "CI/verdict amendment only; source counts and artifacts unchanged",
        }
    )
    shutil.copy2(implementation, campaign.dir / "source_gateb_gb2_amendment.py")
    shutil.copy2(bootstrap_source, campaign.dir / "source_bootstrap.py")
    for name, path in source_files.items():
        shutil.copy2(path, campaign.dir / f"source_{name}")

    verdict = corrected_verdict(
        beam,
        bposd,
        draws=draws,
        seed=seed,
        band=band,
    )
    result = {
        "revision": "GB2",
        "source_campaign": source_campaign.name,
        "raw_counts": {
            "beam8": {"failures": int(beam["failures"]), "shots": int(beam["shots"])},
            "bposd": {"failures": int(bposd["failures"]), "shots": int(bposd["shots"])},
        },
        "source_verdict_superseded": summary.get("verdict"),
        "corrected_verdict": verdict,
    }
    campaign.append_result(result)
    campaign.close(result)
    return campaign.dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=None)
    args = parser.parse_args()
    print(amend(args.campaign, root=args.root))


if __name__ == "__main__":
    main()
