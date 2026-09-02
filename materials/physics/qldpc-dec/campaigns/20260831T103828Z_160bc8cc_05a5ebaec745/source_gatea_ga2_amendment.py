"""Combine a corrected Gate-A accuracy campaign with its GA2 X-row replacement.

The six-row source campaign remains immutable. Its p=3e-4 row labeled X used
an old loader that silently fell back to the derived Z-memory circuit. This
amendment selects the focused, basis-correct replacement for that one row and
reuses the other five corrected-prior rows without decoding again.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from pathlib import Path

from qldpc_dec.bootstrap import wilson_interval
from qldpc_dec.runner import CAMPAIGNS_ROOT, Campaign
from qldpc_dec.seeds import sampling_seed

ROUNDS = 12
BASE_SEED = 20260829
SHOTS = 100_000
BAD_X_CIRCUIT_SHA256 = "c179efd32e7244460c4b59d0a4d200e753469214efc692c389cb3bb4b2523439"
GOOD_X_CIRCUIT_SHA256 = "300f8625ac81c9bf35f16b85b3bebf6aa087947650209387ab0f688168f372a2"
GOOD_X_DEM_SHA256 = "d1d6fdb838d407155ec1a5940b8e1d32b7f2c32d5661dea84dad3feb8f4d4579"
BPOSD_PRIOR_MODEL = "merged DEM heterogeneous error_channel"
CORRECTED_BP_SOURCE_SHA256 = "3a753f46bf404b4bcf07111bfe22c76bf4f3c901b758caa65738e6c6b28055aa"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _load_json_object(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise TypeError(f"expected JSON object in {path}")
    return value


def _load_closed_campaign(path: Path) -> tuple[dict, dict, dict, list[dict]]:
    required = ("manifest.json", "summary.json", "inventory.json", "results.json.gz")
    missing = [name for name in required if not (path / name).is_file()]
    if missing:
        raise FileNotFoundError(f"campaign is not closed; missing {missing}: {path}")
    manifest = _load_json_object(path / "manifest.json")
    summary = _load_json_object(path / "summary.json")
    inventory = _load_json_object(path / "inventory.json")
    if inventory.get("policy") != "immutable-after-close":
        raise ValueError(f"campaign lacks immutable-after-close policy: {path}")
    if set(inventory["files"]) != {item.name for item in path.iterdir()}:
        raise ValueError(f"campaign inventory differs from directory contents: {path}")
    with gzip.open(path / "results.json.gz", "rt") as stream:
        results = json.load(stream)
    if not isinstance(results, list) or len(results) != inventory.get("results_count"):
        raise ValueError(f"campaign result count differs from inventory: {path}")
    source_hashes = manifest.get("source_sha256", {})
    if not isinstance(source_hashes, dict):
        raise TypeError(f"source_sha256 is not an object: {path}")
    for name, expected in source_hashes.items():
        if sha256_file(path / name) != expected:
            raise ValueError(f"source snapshot hash mismatch: {path / name}")
    return manifest, summary, inventory, results


def _validate_row(row: dict, *, p: float, basis: str, shots: int) -> None:
    if (row.get("p"), row.get("basis"), row.get("shots")) != (p, basis, shots):
        raise ValueError(f"unexpected row identity: {row}")
    failures = int(row["failures"])
    raw = failures / shots
    raw_ci = list(wilson_interval(failures, shots))
    if row.get("raw_logical_failure_rate") != raw or row.get("raw_ler_ci95") != raw_ci:
        raise ValueError(f"raw-rate fields are inconsistent: p={p} basis={basis}")
    if row.get("ler") != raw / ROUNDS:
        raise ValueError(f"per-round LER is inconsistent: p={p} basis={basis}")
    if row.get("ler_ci95") != [value / ROUNDS for value in raw_ci]:
        raise ValueError(f"per-round interval is inconsistent: p={p} basis={basis}")
    if row.get("sampling_seed") != sampling_seed(BASE_SEED, p, basis, "bp30+osd"):
        raise ValueError(f"sampling seed is inconsistent: p={p} basis={basis}")
    if row.get("bposd_prior_model") != BPOSD_PRIOR_MODEL:
        raise ValueError(f"row lacks the corrected BP prior model: p={p} basis={basis}")


def amend(primary: Path, replacement: Path, root: Path | None = None) -> Path:
    primary = primary.resolve()
    replacement = replacement.resolve()
    primary_manifest, primary_summary, primary_inventory, primary_rows = (
        _load_closed_campaign(primary)
    )
    replacement_manifest, replacement_summary, replacement_inventory, replacement_rows = (
        _load_closed_campaign(replacement)
    )

    primary_config = primary_manifest.get("config")
    replacement_config = replacement_manifest.get("config")
    if not isinstance(primary_config, dict) or not isinstance(replacement_config, dict):
        raise TypeError("source campaign config is not an object")
    if primary_config != {
        "run": "gateA_accuracy",
        "shots": SHOTS,
        "base_seed": BASE_SEED,
        "bposd_prior_model": BPOSD_PRIOR_MODEL,
    }:
        raise ValueError(f"unexpected primary campaign config: {primary_config}")
    expected_replacement = {
        "run": "gateA_accuracy_GA2_X_p3e-4_replacement",
        "p": 3e-4,
        "basis": "X",
        "shots": SHOTS,
        "base_seed": BASE_SEED,
        "bposd_prior_model": BPOSD_PRIOR_MODEL,
        "replaces": f"{primary.name} p=3e-4 X row",
    }
    if replacement_config != expected_replacement:
        raise ValueError(f"unexpected replacement campaign config: {replacement_config}")
    for label, manifest in (
        ("primary", primary_manifest),
        ("replacement", replacement_manifest),
    ):
        source_hashes = manifest.get("source_sha256")
        if (
            not isinstance(source_hashes, dict)
            or source_hashes.get("source_qldpc_dec_bp_osd.py")
            != CORRECTED_BP_SOURCE_SHA256
        ):
            raise ValueError(f"{label} campaign lacks the audited corrected BP source")
    if len(primary_rows) != 6 or len(replacement_rows) != 1:
        raise ValueError("expected six primary rows and one replacement row")

    by_point: dict[tuple[float, str], dict] = {}
    for row in primary_rows:
        key = (float(row["p"]), str(row["basis"]))
        if key in by_point:
            raise ValueError(f"duplicate primary row {key}")
        by_point[key] = row
    expected_keys = {(p, basis) for p in (3e-4, 5e-4, 1e-3) for basis in ("X", "Z")}
    if set(by_point) != expected_keys:
        raise ValueError(f"primary points differ from the frozen grid: {set(by_point)}")
    for (p, basis), row in by_point.items():
        _validate_row(row, p=p, basis=basis, shots=SHOTS)

    invalid_x = by_point[(3e-4, "X")]
    valid_z = by_point[(3e-4, "Z")]
    if invalid_x["circuit_sha256"] != BAD_X_CIRCUIT_SHA256:
        raise ValueError("primary p=3e-4 X row does not carry the known wrong Z circuit")
    if invalid_x["circuit_sha256"] != valid_z["circuit_sha256"]:
        raise ValueError("primary p=3e-4 X and Z rows do not expose the documented fallback")

    corrected_x = replacement_rows[0]
    _validate_row(corrected_x, p=3e-4, basis="X", shots=SHOTS)
    if corrected_x["circuit_sha256"] != GOOD_X_CIRCUIT_SHA256:
        raise ValueError("replacement does not use the basis-correct X circuit")
    if corrected_x["dem_sha256"] != GOOD_X_DEM_SHA256:
        raise ValueError("replacement does not use the basis-correct X DEM")
    if corrected_x["sampling_seed"] != invalid_x["sampling_seed"]:
        raise ValueError("replacement did not preserve the original X sampling seed")

    selected: list[dict] = []
    for p in (3e-4, 5e-4, 1e-3):
        for basis in ("X", "Z"):
            use_replacement = (p, basis) == (3e-4, "X")
            source = replacement if use_replacement else primary
            row = corrected_x if use_replacement else by_point[(p, basis)]
            selected_row = dict(row)
            selected_row["source_campaign"] = source.name
            selected_row["selection"] = (
                "GA2 basis-correct replacement" if use_replacement else "unchanged valid primary row"
            )
            selected.append(selected_row)

    input_files: dict[str, Path] = {}
    for prefix, source in (("primary", primary), ("replacement", replacement)):
        for name in ("manifest.json", "summary.json", "inventory.json", "results.json.gz"):
            input_files[f"{prefix}_{name}"] = source / name
    input_sha256 = {name: sha256_file(path) for name, path in input_files.items()}
    config = {
        "run": "gateA_accuracy_GA2_amendment",
        "revision": "pre_statement.md GA1-GA2",
        "primary_campaign": primary.name,
        "replacement_campaign": replacement.name,
        "shots_per_point": SHOTS,
        "base_seed": BASE_SEED,
        "decoder_rerun": False,
        "selected_rows": 6,
    }
    campaign = Campaign(config, root=root or CAMPAIGNS_ROOT())
    implementation = Path(__file__).resolve()
    qldpc_source = implementation.parent / "qldpc_dec"
    source_files = {
        "source_gatea_ga2_amendment.py": implementation,
        "source_qldpc_dec_run_gate.py": qldpc_source / "run_gate.py",
        "source_qldpc_dec_bp_osd.py": qldpc_source / "bp_osd.py",
        "source_qldpc_dec_circuits.py": qldpc_source / "circuits.py",
        "source_qldpc_dec_seeds.py": qldpc_source / "seeds.py",
        "source_qldpc_dec_runner.py": qldpc_source / "runner.py",
        "pre_statement.md": implementation.parent.parent / "pre_statement.md",
    }
    for name, source in source_files.items():
        shutil.copy2(source, campaign.dir / name)
    for name, source in input_files.items():
        shutil.copy2(source, campaign.dir / name)
    campaign.write_manifest(
        {
            "primary_closed_utc": primary_inventory["closed_utc"],
            "replacement_closed_utc": replacement_inventory["closed_utc"],
            "input_sha256": input_sha256,
            "source_sha256": {
                name: sha256_file(campaign.dir / name) for name in source_files
            },
            "invalidated_primary_summary": primary_summary,
            "replacement_summary": replacement_summary,
        }
    )
    for row in selected:
        campaign.append_result(row)

    totals: dict[str, dict[str, float | int]] = {}
    for p in (3e-4, 5e-4, 1e-3):
        x = next(row for row in selected if row["p"] == p and row["basis"] == "X")
        z = next(row for row in selected if row["p"] == p and row["basis"] == "Z")
        totals[f"{p:g}"] = {
            "x_failures": int(x["failures"]),
            "z_failures": int(z["failures"]),
            "total_raw_ler_x_plus_z": float(x["raw_logical_failure_rate"])
            + float(z["raw_logical_failure_rate"]),
            "total_per_round_ler_x_plus_z": float(x["ler"]) + float(z["ler"]),
        }
    summary = {
        "gate": "A-accuracy-companion-GA2-amended",
        "status": "DESCRIPTIVE_NO_PUBLISHED_ACCURACY_TARGET",
        "selected_rows": selected,
        "total_ler_by_p": totals,
        "invalidated_primary_row": {
            "p": 3e-4,
            "basis_label": "X",
            "actual_circuit_sha256": invalid_x["circuit_sha256"],
            "reason": "old loader silently fell back to the derived Z-memory circuit",
        },
    }
    campaign.close(summary)
    return campaign.dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--replacement", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=None)
    args = parser.parse_args()
    print(amend(args.primary, args.replacement, root=args.root))


if __name__ == "__main__":
    main()
