"""Repair the omitted p=5e-4 verdict in a closed Gate-A timing campaign.

Revision GA3 changes no timing sample. It validates and copies the immutable
source campaign, then recomputes both published-reference ratios using the
canonical timing-key helper in ``qldpc_dec.run_gate``.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from pathlib import Path

from qldpc_dec.run_gate import (
    BPOSD_PRIOR_MODEL,
    GATE_A_TIMING_REFERENCE,
    gate_a_timing_verdict,
)
from qldpc_dec.runner import CAMPAIGNS_ROOT, Campaign
from qldpc_dec.seeds import sampling_seed

BASE_SEED = 20260829
NUM_TIMED = 10_000
SOURCE_RUN_GATE_SHA256 = "50521a69e378b1847d0dcd8c57884141cef8eb6dc84d052c4f0fa872bda96747"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _load_object(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise TypeError(f"expected JSON object in {path}")
    return value


def amend(source: Path, root: Path | None = None) -> Path:
    source = source.resolve()
    source_files = {
        name: source / name
        for name in ("manifest.json", "summary.json", "inventory.json", "results.json.gz")
    }
    missing = [name for name, path in source_files.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"source timing campaign is not closed; missing {missing}")
    manifest = _load_object(source_files["manifest.json"])
    summary = _load_object(source_files["summary.json"])
    inventory = _load_object(source_files["inventory.json"])
    if inventory.get("policy") != "immutable-after-close":
        raise ValueError("source timing campaign lacks immutable-after-close policy")
    if set(inventory["files"]) != {item.name for item in source.iterdir()}:
        raise ValueError("source timing campaign inventory differs from directory contents")
    with gzip.open(source_files["results.json.gz"], "rt") as stream:
        rows = json.load(stream)
    if not isinstance(rows, list) or len(rows) != 2 or inventory.get("results_count") != 2:
        raise ValueError("source timing campaign must contain exactly two rows")

    config = manifest.get("config")
    expected_config = {
        "run": "gateA_timing",
        "num_timed": NUM_TIMED,
        "base_seed": BASE_SEED,
        "bposd_prior_model": BPOSD_PRIOR_MODEL,
    }
    if config != expected_config:
        raise ValueError(f"unexpected source timing config: {config}")
    source_hashes = manifest.get("source_sha256")
    if not isinstance(source_hashes, dict):
        raise TypeError("source timing manifest lacks source hashes")
    for name, expected in source_hashes.items():
        if sha256_file(source / name) != expected:
            raise ValueError(f"source snapshot hash mismatch: {name}")
    if source_hashes.get("source_qldpc_dec_run_gate.py") != SOURCE_RUN_GATE_SHA256:
        raise ValueError("source campaign does not contain the audited GA3-key-bug runner")
    pre_statement_hash = manifest.get("pre_statement_sha256")
    if pre_statement_hash != sha256_file(source / "pre_statement.md"):
        raise ValueError("source pre-statement snapshot hash mismatch")

    by_key: dict[str, dict] = {}
    for row in rows:
        p = float(row["p"])
        key = f"{p:g}"
        if key in by_key:
            raise ValueError(f"duplicate source timing point {key}")
        if row.get("arm") != "bp30+osd-timing" or row.get("basis") != "Z":
            raise ValueError(f"unexpected timing row identity: {row}")
        if row.get("num_timed") != NUM_TIMED:
            raise ValueError(f"unexpected timing sample count at {key}")
        if row.get("sampling_seed") != sampling_seed(BASE_SEED, p, "Z", "timing"):
            raise ValueError(f"unexpected timing sampling seed at {key}")
        if row.get("bposd_prior_model") != BPOSD_PRIOR_MODEL:
            raise ValueError(f"timing row lacks corrected BP prior model at {key}")
        prior_check = row.get("bposd_prior_check")
        if (
            not isinstance(prior_check, dict)
            or prior_check.get("configured_channel_exact") is not True
            or prior_check.get("columns") != 8784
            or prior_check.get("distinct_probabilities") != 9
        ):
            raise ValueError(f"timing row lacks exact heterogeneous-prior audit at {key}")
        by_key[key] = row
    if set(by_key) != set(GATE_A_TIMING_REFERENCE):
        raise ValueError(f"source timing points are incomplete: {sorted(by_key)}")

    measured = {
        key: {
            "mean_ms": float(by_key[key]["mean_ms"]),
            "p999_ms": float(by_key[key]["p999_ms"]),
        }
        for key in GATE_A_TIMING_REFERENCE
    }
    if summary.get("measured") != measured:
        raise ValueError("source summary measurements differ from source timing rows")
    source_reference = summary.get("published_ref")
    source_verdict = summary.get("verdict")
    if (
        not isinstance(source_reference, dict)
        or "5e-04" not in source_reference
        or "0.0005" in source_reference
        or not isinstance(source_verdict, dict)
        or set(source_verdict) != {"0.001"}
    ):
        raise ValueError("source summary does not expose the documented GA3 key mismatch")

    corrected_verdict = gate_a_timing_verdict(measured)
    if not all(bool(point["within_factor2"]) for point in corrected_verdict.values()):
        outcome = "NOT_REPRODUCED_at_least_one_timing_point_outside_factor2"
    else:
        outcome = "REPRODUCED_all_timing_points_within_factor2"

    input_sha256 = {name: sha256_file(path) for name, path in source_files.items()}
    amendment_config = {
        "run": "gateA_timing_GA3_amendment",
        "revision": "pre_statement.md GA3",
        "source_campaign": source.name,
        "source_config_hash12": manifest.get("config_hash12"),
        "decoder_rerun": False,
        "timing_samples_reused": 2 * NUM_TIMED,
    }
    campaign = Campaign(amendment_config, root=root or CAMPAIGNS_ROOT())
    implementation = Path(__file__).resolve()
    qldpc_source = implementation.parent / "qldpc_dec"
    snapshot_sources = {
        "source_gatea_timing_ga3_amendment.py": implementation,
        "source_qldpc_dec_run_gate.py": qldpc_source / "run_gate.py",
        "source_qldpc_dec_runner.py": qldpc_source / "runner.py",
        "pre_statement.md": implementation.parent.parent / "pre_statement.md",
    }
    for name, path in snapshot_sources.items():
        shutil.copy2(path, campaign.dir / name)
    for name, path in source_files.items():
        shutil.copy2(path, campaign.dir / f"source_{name}")
    campaign.write_manifest(
        {
            "source_campaign_closed_utc": inventory["closed_utc"],
            "source_file_sha256": input_sha256,
            "source_sha256": {
                name: sha256_file(campaign.dir / name) for name in snapshot_sources
            },
            "scope": "reference-key verdict repair only; no timing sample rerun",
        }
    )
    result = {
        "revision": "GA3",
        "source_campaign": source.name,
        "measured": measured,
        "published_ref": GATE_A_TIMING_REFERENCE,
        "source_verdict_superseded": source_verdict,
        "corrected_verdict": corrected_verdict,
        "outcome": outcome,
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
