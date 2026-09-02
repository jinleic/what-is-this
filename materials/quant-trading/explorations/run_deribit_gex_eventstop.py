#!/usr/bin/env python3
"""Run one event-count-stopped Deribit taker-GEX fast-screen trial on frozen prefixes.

The parent screen rules (GEX formula, RV completeness, sample minima, gates,
earliest-stop search, truth table) are executed by the hash-pinned parent module
through ``replicate_accelerated_intraday.run_deribit``. This runner adds only the
preregistered window derivation: the window starts at the frozen contract start
and ends at the last 15-minute boundary covered by every required stream inside
the byte prefixes pinned in the contract. There is no wall-clock deadline.

Lifecycle: verify pins -> preflight outputs -> prepare everything in memory ->
re-verify pins -> publish snapshot exclusively -> compute -> re-verify pins ->
publish results exclusively. Nothing is written before every precondition holds,
and no output path is ever overwritten.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import replicate_accelerated_intraday as BASE
import replicate_accelerated_intraday_corrected as CORRECTED


WORKSPACE_ROOT = BASE.WORKSPACE_ROOT
EXPLORATIONS = BASE.EXPLORATIONS
RUNNER_PATH = Path(__file__).resolve()
DEFAULT_CONTRACT = EXPLORATIONS / "deribit-gex-eventstop-contract.json"
BLOCK_MS = 15 * 60 * 1000
SCREEN = "deribit_intraday_taker_gex"
CONTRACT_STATUS = "frozen-before-eventstop-outcome-analysis"
CLASSIFICATION = (
    "event-count-stopped fast-screen trial on append-only public capture; "
    "outcome-blind at registration; not confirmation, not a rescue, no alpha claim"
)
PIN_SECTIONS = ("frozen_rule_files", "baseline_artifact_hashes", "governing_authority")


def floor_block(timestamp_ms: int) -> int:
    return timestamp_ms - timestamp_ms % BLOCK_MS


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _resolve(record: dict) -> Path:
    return BASE.resolve_workspace_path(record["path"])


def stream_times(records: dict[str, list[dict]]) -> dict[str, list[int]]:
    """Authoritative exchange timestamps per required Deribit stream."""
    times: dict[str, list[int]] = {
        "deribit_option_trade": [],
        "btc_usd_index": [],
        "eth_usd_index": [],
    }
    for path, rows in records.items():
        name = Path(path).name
        if name == "trades.jsonl":
            for row in rows:
                timestamp = (row.get("trade") or {}).get("timestamp")
                if _is_int(timestamp):
                    times["deribit_option_trade"].append(timestamp)
        elif name == "index.jsonl":
            for row in rows:
                data = row.get("data") or {}
                key = f"{data.get('index_name')}_index"
                timestamp = data.get("timestamp")
                if key in times and _is_int(timestamp):
                    times[key].append(timestamp)
    return times


def derive_window_end(records: dict[str, list[dict]], start_ms: int) -> dict[str, Any]:
    """end = last block boundary <= min(latest exchange ts) over required streams.

    Uses timestamps only; fails closed on any missing stream. Never reads a
    price, GEX, RV, or gate value.
    """
    times = stream_times(records)
    missing = [name for name, values in times.items() if not values]
    if missing:
        raise ValueError(f"required stream missing from frozen prefixes: {missing}")
    latest = {name: max(values) for name, values in times.items()}
    common_cutoff_ms = min(latest.values())
    end_ms = floor_block(common_cutoff_ms)
    if end_ms <= start_ms:
        raise ValueError(
            f"frozen prefixes end before the first complete block: cutoff {common_cutoff_ms} <= start {start_ms}"
        )
    streams = {name: CORRECTED._coverage(values, end_ms) for name, values in times.items()}
    failed = [name for name, report in streams.items() if not report["covers_end"]]
    if failed:
        raise ValueError(f"required stream does not cover derived end: {failed}")
    return {
        "rule": "largest 15-minute-aligned boundary <= min over required streams of the latest exchange timestamp inside the frozen prefixes",
        "required_streams": sorted(times),
        "per_stream_latest_ms": latest,
        "per_stream_latest_utc": {name: BASE.iso_ms(value) for name, value in latest.items()},
        "common_cutoff_ms": common_cutoff_ms,
        "common_cutoff_utc": BASE.iso_ms(common_cutoff_ms),
        "end_ms": end_ms,
        "end_utc_exclusive": BASE.iso_ms(end_ms),
        "all_required_streams_cover_end": True,
        "streams": streams,
    }


def validate_feasibility(
    start_ms: int, end_ms: int, required_transitions: int, *, not_before_ms: int
) -> dict[str, Any]:
    """Prove the registered transition minimum is reachable before any outcome access."""
    if start_ms % BLOCK_MS or end_ms % BLOCK_MS:
        raise ValueError("window boundaries must be aligned to parent 15-minute blocks")
    if end_ms <= start_ms:
        raise ValueError("window end must be after start")
    if start_ms < not_before_ms:
        raise ValueError(f"window overlaps a prior evaluated window: {start_ms} < {not_before_ms}")
    window_blocks = (end_ms - start_ms) // BLOCK_MS
    complete_signal_blocks = window_blocks - 1
    max_transitions = 2 * complete_signal_blocks
    report = {
        "block_ms": BLOCK_MS,
        "window_duration_ms": end_ms - start_ms,
        "window_minutes": (end_ms - start_ms) // 60_000,
        "window_blocks": window_blocks,
        "complete_signal_blocks": complete_signal_blocks,
        "max_complete_asset_block_transitions": max_transitions,
        "required_complete_asset_block_transitions": required_transitions,
        "not_before_ms": not_before_ms,
        "not_before_utc": BASE.iso_ms(not_before_ms),
        "nonoverlap_with_prior_window": True,
        "registered_stop_reachable": max_transitions >= required_transitions,
    }
    if not report["registered_stop_reachable"]:
        raise ValueError(
            f"registered sample stop is structurally unreachable: {max_transitions}/{required_transitions} transitions"
        )
    return report


def verify_pins(contract: dict) -> dict[str, dict[str, Any]]:
    """Every pinned file must still hash as frozen. Raises on any drift."""
    runner = contract["runner"]
    if _resolve(runner) != RUNNER_PATH or BASE.sha256_file(RUNNER_PATH) != runner["sha256"]:
        raise ValueError("event-stop runner hash mismatch")
    observed: dict[str, dict[str, Any]] = {}
    for key in ("parent_contract", "corrected_contract"):
        CORRECTED._verify_hash_record(contract["lineage"][key])
    for section in PIN_SECTIONS:
        for key, record in contract[section].items():
            observed_sha256 = BASE.sha256_file(_resolve(record))
            if observed_sha256 != record["sha256"]:
                raise ValueError(
                    f"frozen file changed: {record['path']}: {observed_sha256} != {record['sha256']}"
                )
            observed[f"{section}.{key}"] = {
                "path": record["path"],
                "expected_sha256": record["sha256"],
                "observed_sha256": observed_sha256,
                "unchanged": True,
            }
    return observed


def load_and_verify_contract(path: Path) -> tuple[dict, dict, str, int, int]:
    contract_bytes = path.read_bytes()
    contract_sha256 = BASE.sha256_bytes(contract_bytes)
    contract = json.loads(contract_bytes)
    if contract.get("status") != CONTRACT_STATUS:
        raise ValueError("event-stop contract is not frozen")
    verify_pins(contract)
    parent = json.loads(_resolve(contract["lineage"]["parent_contract"]).read_text(encoding="utf-8"))
    if parent["screens"][SCREEN] != contract["screen_definition"]:
        raise ValueError("screen_definition is not a verbatim copy of the parent screen")
    corrected = json.loads(
        _resolve(contract["lineage"]["corrected_contract"]).read_text(encoding="utf-8")
    )
    start_ms = BASE.parse_iso_ms(contract["observation_window"]["start_utc_inclusive"])
    not_before_ms = BASE.parse_iso_ms(corrected["observation_window"]["end_utc_exclusive"])
    if start_ms != not_before_ms:
        raise ValueError("window must start exactly at the corrected window's exclusive end")
    return contract, parent, contract_sha256, start_ms, not_before_ms


def load_manifest(run_dir: Path, contract_sha256: str) -> dict:
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("prereg_sha256") != contract_sha256:
        raise ValueError("campaign manifest prereg_sha256 does not bind this contract")
    if manifest.get("status") != "RUNNING":
        raise ValueError("campaign is not RUNNING")
    return manifest


def read_contract_prefix(item: dict) -> tuple[dict, list[dict]]:
    path = _resolve(item)
    prefix, _ = BASE.read_exact_prefix(path, item["prefix_bytes"])
    observed = BASE.sha256_bytes(prefix)
    if observed != item["sha256"]:
        raise ValueError(f"{path}: frozen prefix hash mismatch {observed} != {item['sha256']}")
    records, parse = BASE.json_records_from_prefix(prefix, path)
    if path.name == "instruments.jsonl":
        clocks = CORRECTED.instrument_clock_provenance(records)
    else:
        timestamps = [
            value
            for row in records
            if (value := BASE.authoritative_timestamp(path, row)) is not None
        ]
        clocks = {
            "exchange_clock": {
                "field": "trade.timestamp" if path.name == "trades.jsonl" else "data.timestamp",
                "first_ms": min(timestamps) if timestamps else None,
                "first_utc": BASE.iso_ms(min(timestamps)) if timestamps else None,
                "latest_ms": max(timestamps) if timestamps else None,
                "latest_utc": BASE.iso_ms(max(timestamps)) if timestamps else None,
            }
        }
    report = {
        "path": item["path"],
        "prefix_bytes": item["prefix_bytes"],
        "sha256": item["sha256"],
        "parse": parse,
        **clocks,
    }
    return report, records


def prepare(contract: dict, parent: dict, contract_sha256: str, start_ms: int, not_before_ms: int) -> dict[str, Any]:
    """Everything that must hold before any byte is written, computed in memory only."""
    reports: list[dict] = []
    records: dict[str, list[dict]] = {}
    for item in contract["frozen_input_prefixes"]:
        report, parsed = read_contract_prefix(item)
        reports.append(report)
        records[item["path"]] = parsed
    snapshot = {
        "schema_version": 1,
        "created_at_utc": BASE.now_utc(),
        "contract_sha256": contract_sha256,
        "snapshot_rule": "byte prefixes pinned in the contract before campaign init; hash verified; one trailing partial line ignored; malformed complete lines fail",
        "metric_order": "manifest published after every precondition held and before any GEX, RV, correlation, gate, or verdict computation",
        "files": reports,
    }
    window = derive_window_end(records, start_ms)
    feasibility = validate_feasibility(
        start_ms,
        window["end_ms"],
        parent["screens"][SCREEN]["sample_stop"]["minimum_complete_asset_block_transitions"],
        not_before_ms=not_before_ms,
    )
    provenance = CORRECTED.instrument_file_provenance(records, snapshot)
    return {
        "records": records,
        "snapshot": snapshot,
        "window": window,
        "feasibility": feasibility,
        "instrument_provenance": provenance,
    }


def encode_json(value: dict) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=False, ensure_ascii=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def write_bytes_exclusive(path: Path, payload: bytes) -> str:
    """Create-only, fsynced publication; never replaces an existing file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(fd, view)
            view = view[written:]
        os.fsync(fd)
    finally:
        os.close(fd)
    return BASE.sha256_bytes(payload)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    contract_path = args.contract.resolve()
    run_dir = args.run_dir.resolve()

    contract, parent, contract_sha256, start_ms, not_before_ms = load_and_verify_contract(contract_path)
    manifest = load_manifest(run_dir, contract_sha256)
    outputs = contract["outputs"]
    snapshot_path = run_dir / outputs["snapshot_filename"]
    results_path = run_dir / outputs["results_filename"]
    owner_copy = BASE.resolve_workspace_path(outputs["owner_dir_copy"])
    if snapshot_path.exists():
        raise FileExistsError(
            f"existing snapshot {snapshot_path}: a prior invocation reached publication; close the campaign CRASHED"
        )
    CORRECTED._preflight_outputs([results_path, owner_copy])

    prepared = prepare(contract, parent, contract_sha256, start_ms, not_before_ms)
    end_ms = prepared["window"]["end_ms"]

    verify_pins(contract)
    deribit_module = BASE.load_module("eventstop_parent_deribit", BASE.DERIBIT_MODULE)
    snapshot_sha256 = write_bytes_exclusive(snapshot_path, encode_json(prepared["snapshot"]))
    result = BASE.run_deribit(parent, prepared["records"], start_ms, end_ms, deribit_module)
    result["classification"] = CLASSIFICATION
    pins_at_publication = verify_pins(contract)

    common = {
        "eventstop_contract": {
            "path": str(contract_path.relative_to(WORKSPACE_ROOT)),
            "sha256": contract_sha256,
            "created_at_utc": contract["created_at_utc"],
            "classification": contract["classification"],
        },
        "campaign_manifest": manifest,
        "lineage": contract["lineage"],
        "stopping_rule": contract["stopping_rule"],
        "claim_boundary": contract["claim_boundary"],
        "governance": contract["governance"],
        "integrity_preconditions": {
            "window_derivation": prepared["window"],
            "window_feasibility": prepared["feasibility"],
            "instrument_clock_provenance": prepared["instrument_provenance"],
        },
        "input_snapshot": {
            "path": str(snapshot_path.relative_to(WORKSPACE_ROOT)),
            "sha256": snapshot_sha256,
            "file_count": len(prepared["snapshot"]["files"]),
        },
        "pins_verified_at_publication": pins_at_publication,
    }
    payload = encode_json({**common, **result})
    results_sha256 = write_bytes_exclusive(results_path, payload)
    copy_sha256 = write_bytes_exclusive(owner_copy, payload)
    if copy_sha256 != results_sha256:
        raise ValueError("owner copy hash differs from published results")
    print(
        json.dumps(
            {
                "contract_sha256": contract_sha256,
                "run_id": manifest["run_id"],
                "window": {
                    "start_utc_inclusive": BASE.iso_ms(start_ms),
                    "end_utc_exclusive": prepared["window"]["end_utc_exclusive"],
                },
                "sample_stop": result["sample_stop"]["stop_outcome_end_utc"],
                "sample_all_met": result["sample_conditions"]["all_met"],
                "verdict": result["decision"]["verdict"],
                "results": str(results_path),
                "results_sha256": results_sha256,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
