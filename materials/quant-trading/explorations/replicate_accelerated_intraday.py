#!/usr/bin/env python3
"""Run the frozen 2026-08-31 exact-rule chronological replication.

The script snapshots append-only public JSONL inputs before metric computation,
filters the frozen half-open UTC window at ingest, and delegates every metric,
gate, and truth rule to the parent screen implementations whose hashes are
pinned in accelerated-intraday-replication-contract.json.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
QUANT_ROOT = WORKSPACE_ROOT / "quant-trading"
EXPLORATIONS = QUANT_ROOT / "explorations"
RUNNER_PATH = Path(__file__).resolve()
DEFAULT_CONTRACT = EXPLORATIONS / "accelerated-intraday-replication-contract.json"
DEFAULT_SNAPSHOT = (
    QUANT_ROOT
    / "data/derived/accelerated-intraday-replication-20260831/snapshot.json"
)
DEFAULT_COMBINED = EXPLORATIONS / "accelerated-intraday-replication-results.json"
LIQUIDATION_DIR = QUANT_ROOT / "data/raw/focused/live/liquidations-persistent"
DERIBIT_DIR = QUANT_ROOT / "data/raw/focused/live/deribit"
LIQUIDATION_MODULE = EXPLORATIONS / "focused-liquidation-overlay/fast_intraday.py"
DERIBIT_MODULE = EXPLORATIONS / "focused-deribit-gex/fast_intraday.py"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, limit: int | None = None) -> str:
    digest = hashlib.sha256()
    remaining = limit
    with path.open("rb") as handle:
        while True:
            chunk_size = 1024 * 1024 if remaining is None else min(1024 * 1024, remaining)
            if chunk_size <= 0:
                break
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
            if remaining is not None:
                remaining -= len(chunk)
    return digest.hexdigest()


def parse_iso_ms(text: str) -> int:
    parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError(f"timestamp lacks timezone: {text}")
    return int(round(parsed.timestamp() * 1000))


def iso_ms(timestamp_ms: int | None) -> str | None:
    if timestamp_ms is None:
        return None
    parsed = datetime.fromtimestamp(timestamp_ms / 1000, timezone.utc)
    return parsed.strftime("%Y-%m-%dT%H:%M:%S.") + f"{timestamp_ms % 1000:03d}Z"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def in_half_open(timestamp_ms: int, start_ms: int, end_ms: int) -> bool:
    return start_ms <= timestamp_ms < end_ms


def resolve_workspace_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE_ROOT / path


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def json_records_from_prefix(prefix: bytes, path: Path) -> tuple[list[dict], dict]:
    trailing_partial_bytes = 0
    complete = prefix
    if prefix and not prefix.endswith(b"\n"):
        split_at = prefix.rfind(b"\n")
        if split_at < 0:
            trailing_partial_bytes = len(prefix)
            complete = b""
        else:
            trailing_partial_bytes = len(prefix) - split_at - 1
            complete = prefix[: split_at + 1]

    records: list[dict] = []
    blank_lines = 0
    malformed: list[dict] = []
    for line_number, raw in enumerate(complete.splitlines(), 1):
        if not raw.strip():
            blank_lines += 1
            continue
        try:
            row = json.loads(raw)
        except json.JSONDecodeError as exc:
            malformed.append({"line": line_number, "error": str(exc)})
            continue
        if not isinstance(row, dict):
            malformed.append({"line": line_number, "error": "JSON value is not an object"})
            continue
        records.append(row)
    if malformed:
        raise ValueError(f"{path}: malformed complete JSON lines: {malformed[:3]}")
    return records, {
        "complete_lines": len(complete.splitlines()),
        "blank_lines": blank_lines,
        "records": len(records),
        "trailing_partial_line_bytes_ignored": trailing_partial_bytes,
        "malformed_complete_lines": 0,
    }


def authoritative_timestamp(path: Path, row: dict) -> int | None:
    if path.name.startswith("liquidation-capture-"):
        if row.get("kind") == "event":
            payload = row.get("payload")
            if isinstance(payload, dict) and isinstance(payload.get("data"), dict):
                payload = payload["data"]
            value = payload.get("E") if isinstance(payload, dict) else None
            return value if isinstance(value, int) and not isinstance(value, bool) else None
        value = row.get("ts_utc")
        if isinstance(value, str):
            try:
                return parse_iso_ms(value)
            except ValueError:
                return None
        return None
    if path.name == "trades.jsonl":
        value = (row.get("trade") or {}).get("timestamp")
    elif path.name == "index.jsonl":
        value = (row.get("data") or {}).get("timestamp")
    elif path.name == "instruments.jsonl":
        value = row.get("received_ts_ms")
    else:
        value = None
    return value if isinstance(value, int) and not isinstance(value, bool) else None


def read_exact_prefix(path: Path, prefix_bytes: int) -> tuple[bytes, os.stat_result]:
    stat = path.stat()
    if stat.st_size < prefix_bytes:
        raise ValueError(f"{path}: current size {stat.st_size} is below frozen prefix {prefix_bytes}")
    with path.open("rb") as handle:
        prefix = handle.read(prefix_bytes)
    if len(prefix) != prefix_bytes:
        raise ValueError(f"{path}: short prefix read {len(prefix)} != {prefix_bytes}")
    return prefix, stat


def snapshot_file(path: Path) -> tuple[dict, list[dict]]:
    stat_before = path.stat()
    prefix, _ = read_exact_prefix(path, stat_before.st_size)
    records, parse = json_records_from_prefix(prefix, path)
    timestamps = [
        timestamp
        for row in records
        if (timestamp := authoritative_timestamp(path, row)) is not None
    ]
    report = {
        "path": str(path.relative_to(WORKSPACE_ROOT)),
        "prefix_bytes": stat_before.st_size,
        "sha256": sha256_bytes(prefix),
        "mtime_ns_at_stat": stat_before.st_mtime_ns,
        "mtime_utc_at_stat": datetime.fromtimestamp(
            stat_before.st_mtime_ns / 1_000_000_000, timezone.utc
        ).isoformat(timespec="microseconds").replace("+00:00", "Z"),
        "snapshot_wall_clock_utc": now_utc(),
        "first_authoritative_timestamp_ms": min(timestamps) if timestamps else None,
        "latest_authoritative_timestamp_ms": max(timestamps) if timestamps else None,
        "first_authoritative_timestamp_utc": iso_ms(min(timestamps)) if timestamps else None,
        "latest_authoritative_timestamp_utc": iso_ms(max(timestamps)) if timestamps else None,
        "parse": parse,
    }
    return report, records


def read_frozen_snapshot(report: dict) -> list[dict]:
    path = resolve_workspace_path(report["path"])
    prefix, _ = read_exact_prefix(path, report["prefix_bytes"])
    observed = sha256_bytes(prefix)
    if observed != report["sha256"]:
        raise ValueError(f"{path}: frozen prefix hash mismatch {observed} != {report['sha256']}")
    records, parse = json_records_from_prefix(prefix, path)
    if parse != report["parse"]:
        raise ValueError(f"{path}: parse summary changed for frozen prefix")
    return records


def collect_input_paths() -> list[Path]:
    paths: list[Path] = []
    for day in ("20260830", "20260831"):
        event_path = LIQUIDATION_DIR / f"liquidation-capture-{day}.jsonl"
        if event_path.is_file():
            paths.append(event_path)
        day_dir = DERIBIT_DIR / f"day{day}"
        for name in ("trades.jsonl", "instruments.jsonl", "index.jsonl"):
            candidate = day_dir / name
            if candidate.is_file():
                paths.append(candidate)
    required = {
        "liquidation": any(p.name.startswith("liquidation-capture-") for p in paths),
        "trades": any(p.name == "trades.jsonl" for p in paths),
        "instruments": any(p.name == "instruments.jsonl" for p in paths),
        "index": any(p.name == "index.jsonl" for p in paths),
    }
    missing = [name for name, present in required.items() if not present]
    if missing:
        raise ValueError(f"required public inputs missing: {missing}")
    return sorted(paths)


def write_json_new(path: Path, value: dict) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite frozen artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=False, ensure_ascii=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def load_or_create_snapshot(contract_path: Path, contract_sha256: str, output: Path) -> tuple[dict, dict[str, list[dict]]]:
    runner_record = {
        "path": str(RUNNER_PATH.relative_to(WORKSPACE_ROOT)),
        "sha256": sha256_file(RUNNER_PATH),
    }
    if output.exists():
        manifest = json.loads(output.read_text(encoding="utf-8"))
        if manifest.get("replication_contract_sha256") != contract_sha256:
            raise ValueError("existing snapshot belongs to a different replication contract")
        if manifest.get("replication_runner") != runner_record:
            raise ValueError("replication runner changed after the snapshot was frozen")
        records = {
            item["path"]: read_frozen_snapshot(item)
            for item in manifest["files"]
        }
        return manifest, records

    reports: list[dict] = []
    records: dict[str, list[dict]] = {}
    for path in collect_input_paths():
        report, parsed = snapshot_file(path)
        reports.append(report)
        records[report["path"]] = parsed
    manifest = {
        "schema_version": 1,
        "created_at_utc": now_utc(),
        "replication_contract": str(contract_path.relative_to(WORKSPACE_ROOT)),
        "replication_contract_sha256": contract_sha256,
        "replication_runner": runner_record,
        "snapshot_rule": "stat size first; hash and parse exactly that prefix; ignore one trailing partial line; fail on malformed complete lines",
        "metric_order": "manifest written before any return, GEX, RV, correlation, gate, or verdict computation",
        "files": reports,
    }
    write_json_new(output, manifest)
    return manifest, records


def verify_hash_record(record: dict) -> None:
    path = resolve_workspace_path(record["path"])
    observed = sha256_file(path)
    if observed != record["sha256"]:
        raise ValueError(f"frozen file changed: {path}: {observed} != {record['sha256']}")


def load_and_verify_contract(path: Path) -> tuple[dict, dict, str, int, int]:
    contract_bytes = path.read_bytes()
    contract_sha256 = sha256_bytes(contract_bytes)
    contract = json.loads(contract_bytes)
    if contract.get("status") != "frozen-before-replication-outcome-analysis":
        raise ValueError("replication contract is not frozen")
    parent_record = contract["parent_contract"]
    parent_path = resolve_workspace_path(parent_record["path"])
    parent_bytes = parent_path.read_bytes()
    if sha256_bytes(parent_bytes) != parent_record["sha256"]:
        raise ValueError("parent contract hash mismatch")
    parent = json.loads(parent_bytes)
    for record in contract["frozen_rule_files"].values():
        verify_hash_record(record)
    for record in contract["baseline_artifact_hashes"].values():
        verify_hash_record(record)
    start_ms = parse_iso_ms(contract["observation_window"]["start_utc_inclusive"])
    end_ms = parse_iso_ms(contract["observation_window"]["end_utc_exclusive"])
    if end_ms - start_ms != 90 * 60 * 1000:
        raise ValueError("replication window is not exactly 90 minutes")
    return contract, parent, contract_sha256, start_ms, end_ms


def force_order_event_ms(row: dict) -> int | None:
    payload = row.get("payload")
    if isinstance(payload, dict) and isinstance(payload.get("data"), dict):
        payload = payload["data"]
    value = payload.get("E") if isinstance(payload, dict) else None
    return value if isinstance(value, int) and not isinstance(value, bool) else None


def filter_liquidation_context(rows: list[dict], end_ms: int) -> list[dict]:
    filtered: list[dict] = []
    for row in rows:
        if row.get("kind") == "event":
            timestamp = force_order_event_ms(row)
        else:
            value = row.get("ts_utc")
            try:
                timestamp = parse_iso_ms(value) if isinstance(value, str) else None
            except ValueError:
                timestamp = None
        if timestamp is not None and timestamp < end_ms:
            filtered.append(row)
    return filtered


def filter_liquidation_onsets(onsets: list[dict], start_ms: int, end_ms: int) -> list[dict]:
    return [
        onset
        for onset in onsets
        if in_half_open(onset["onset_ms"], start_ms, end_ms)
    ]


def filter_liquidation_index(rows: list[dict], start_ms: int, end_ms: int) -> tuple[list[dict], dict]:
    lower = start_ms - 5_000
    kept: list[dict] = []
    counts = {"before_anchor_context": 0, "at_or_after_end": 0, "missing_timestamp": 0}
    for row in rows:
        timestamp = (row.get("data") or {}).get("timestamp")
        if not isinstance(timestamp, int) or isinstance(timestamp, bool):
            counts["missing_timestamp"] += 1
        elif timestamp < lower:
            counts["before_anchor_context"] += 1
        elif timestamp >= end_ms:
            counts["at_or_after_end"] += 1
        else:
            kept.append(row)
    return kept, counts


def filter_deribit_records(records: dict[str, list[dict]], start_ms: int, end_ms: int) -> tuple[dict[str, list[dict]], dict]:
    kept: dict[str, list[dict]] = {"trades": [], "index": [], "instruments": []}
    counts = {
        "excluded_before_start_or_anchor": {"trades": 0, "index": 0, "instruments": 0},
        "excluded_at_or_after_end": {"trades": 0, "index": 0, "instruments": 0},
        "excluded_missing_timestamp": {"trades": 0, "index": 0, "instruments": 0},
    }
    lower_by_kind = {"trades": start_ms, "index": start_ms - 5_000, "instruments": -math.inf}
    for kind in ("trades", "index", "instruments"):
        for row in records[kind]:
            if kind == "trades":
                timestamp = (row.get("trade") or {}).get("timestamp")
            elif kind == "index":
                timestamp = (row.get("data") or {}).get("timestamp")
            else:
                timestamp = row.get("received_ts_ms")
            if not isinstance(timestamp, int) or isinstance(timestamp, bool):
                counts["excluded_missing_timestamp"][kind] += 1
            elif timestamp < lower_by_kind[kind]:
                counts["excluded_before_start_or_anchor"][kind] += 1
            elif timestamp >= end_ms:
                counts["excluded_at_or_after_end"][kind] += 1
            else:
                kept[kind].append(row)
    return kept, counts


def records_by_filename(records: dict[str, list[dict]], filename: str) -> list[dict]:
    combined: list[dict] = []
    for path, rows in sorted(records.items()):
        if Path(path).name == filename:
            for sequence, row in enumerate(rows):
                copy = dict(row)
                copy["_seq"] = sequence
                combined.append(copy)
    return combined


def run_liquidation(parent: dict, records: dict[str, list[dict]], start_ms: int, end_ms: int, module: ModuleType) -> dict:
    event_rows = [
        row
        for path, rows in sorted(records.items())
        if Path(path).name.startswith("liquidation-capture-")
        for row in rows
    ]
    context_rows = filter_liquidation_context(event_rows, end_ms)
    index_rows, index_exclusions = filter_liquidation_index(
        records_by_filename(records, "index.jsonl"), start_ms, end_ms
    )
    cutoff_ms = end_ms - 1
    onset_stats = module.build_onsets(context_rows, cutoff_ms)
    eligible = filter_liquidation_onsets(
        onset_stats["eligible_onsets"], start_ms, end_ms
    )
    index = module.load_index_series(index_rows, cutoff_ms)
    outcomes = module.compute_onset_outcomes(eligible, index, cutoff_ms)
    screen = parent["screens"]["liquidation_btc_eth_post_onset"]
    stop = module.find_earliest_stop(outcomes, screen["sample_stop"], end_ms)
    cohort = stop["cohort"]
    conditions = stop["conditions"]
    gates = module.evaluate_gates(cohort, screen["gate"], conditions)
    truth = module.resolve_truth(stop["found"], gates)
    return {
        "schema_version": 1,
        "screen": "liquidation_btc_eth_post_onset",
        "classification": "outcome-blind retrospective chronological replication",
        "window": {
            "start_utc_inclusive": iso_ms(start_ms),
            "end_utc_exclusive": iso_ms(end_ms),
            "event_rows_before_window_filter": len(event_rows),
            "event_context_rows_before_end": len(context_rows),
            "eligible_onsets_all_context": len(onset_stats["eligible_onsets"]),
            "eligible_onsets_in_window": len(eligible),
            "index_rows_in_anchor_or_window": len(index_rows),
            "index_exclusions": index_exclusions,
        },
        "parent_rules": {
            "sample_stop": screen["sample_stop"],
            "gates": screen["gate"],
            "truth_table": screen["truth_table"],
        },
        "sample_stop": {
            "found": stop["found"],
            "stop_ms": stop["stop_ms"],
            "stop_utc": iso_ms(stop["stop_ms"]),
            "candidate_attempts": stop["tried"],
        },
        "sample_conditions": {
            "all_met": stop["found"],
            "conditions": conditions["conditions"],
            "by_symbol": conditions["by_symbol"],
            "by_side": conditions["by_side"],
            "nonempty_15m_blocks": conditions["nonempty_15m_blocks"],
        },
        "cohort": {
            "size": len(cohort),
            "event_ids": [row["event_id"] for row in cohort],
            "detail": cohort,
        },
        "gates": gates,
        "truth": truth,
        "decision_boundary": "A replication pass cannot reverse the parent falsification or promote alpha.",
    }


def run_deribit(parent: dict, records: dict[str, list[dict]], start_ms: int, end_ms: int, module: ModuleType) -> dict:
    module.RUN = module.load_run_module()
    raw_by_kind = {
        "trades": records_by_filename(records, "trades.jsonl"),
        "instruments": records_by_filename(records, "instruments.jsonl"),
        "index": records_by_filename(records, "index.jsonl"),
    }
    windowed, exclusions = filter_deribit_records(raw_by_kind, start_ms, end_ms)
    all_records = windowed["trades"] + windowed["instruments"] + windowed["index"]
    module.replay_order(all_records)
    deduped, duplicate_drops = module.dedupe(all_records)
    by_kind: dict[str, list[dict]] = {"trades": [], "instruments": [], "index": []}
    for row in deduped:
        record_type = row.get("record_type")
        if record_type == "option_trade":
            by_kind["trades"].append(row)
        elif record_type == "instrument":
            by_kind["instruments"].append(row)
        elif record_type == "index_price":
            by_kind["index"].append(row)

    instruments, metadata_rows = module.build_instrument_map(by_kind)
    exclusion_reasons: dict[str, int] = defaultdict(int)
    gex_rows = []
    for row in by_kind["trades"]:
        trade = row.get("trade") or {}
        computed = module.compute_trade_gex(
            trade, instruments.get(trade.get("instrument_name"))
        )
        gex_rows.append(computed)
        if not computed["ok"]:
            exclusion_reasons[computed["reason"]] += 1
    series, index_counts = module.build_index_series(by_kind["index"])
    aggregate = module.compute_blocks(gex_rows, ts_max=None)
    transitions_all, auxiliary_all = module.build_transitions(aggregate, series)
    transitions_window = [
        row
        for row in transitions_all
        if row["signal_block"] * module.BLOCK_MS >= start_ms
        and row["outcome_end_ms"] <= end_ms
    ]
    screen = parent["screens"]["deribit_intraday_taker_gex"]
    stop_ms, trail = module.find_earliest_stop(
        transitions_window, screen["sample_stop"], end_ms
    )
    transitions = (
        [row for row in transitions_window if row["outcome_end_ms"] <= stop_ms]
        if stop_ms is not None
        else transitions_window
    )
    conditions = module.conditions_at(transitions, screen["sample_stop"])
    gates, extras = module.evaluate_gates(transitions, screen["gate"])
    verdict, reasons = module.resolve_verdict(conditions["conditions_detail"], gates)
    return {
        "schema_version": 1,
        "screen": "deribit_intraday_taker_gex",
        "classification": "outcome-blind retrospective chronological replication",
        "window": {
            "start_utc_inclusive": iso_ms(start_ms),
            "end_utc_exclusive": iso_ms(end_ms),
            "raw_records_by_kind": {name: len(rows) for name, rows in raw_by_kind.items()},
            "filtered_records_by_kind": {name: len(rows) for name, rows in windowed.items()},
            "filter_exclusions": exclusions,
            "records_after_dedup": {name: len(rows) for name, rows in by_kind.items()},
            "duplicate_drops": duplicate_drops,
        },
        "parent_rules": {
            "sample_stop": screen["sample_stop"],
            "gates": screen["gate"],
            "truth_table": screen["truth_table"],
        },
        "reconciliation": {
            "instrument_metadata_rows": metadata_rows,
            "instrument_map_size": len(instruments),
            "gex_rows": len(gex_rows),
            "gex_computable": sum(1 for row in gex_rows if row["ok"]),
            "gex_excluded_by_reason": dict(sorted(exclusion_reasons.items())),
            "index_counts": index_counts,
        },
        "sample_stop": {
            "found": stop_ms is not None,
            "stop_outcome_end_ms": stop_ms,
            "stop_outcome_end_utc": iso_ms(stop_ms),
            "candidate_trail": trail,
        },
        "transitions": {
            "eligible_before_window_end": len(transitions_window),
            "included_count": len(transitions),
            "included": transitions,
            "rv_journal": auxiliary_all["rv_journal"],
            "zero_gex_blocks": auxiliary_all["zero_gex_blocks"],
        },
        "sample_conditions": {
            "all_met": conditions["met"],
            "conditions": conditions["conditions_detail"],
        },
        "gates": {
            "evaluated": gates,
            "pooled_observations": extras["pooled_n"],
            "per_asset_detail": extras["per_asset_spearman_detail"],
        },
        "decision": {
            "verdict": verdict,
            "reasons": reasons,
            "no_returns": True,
            "decision_boundary": "A replication pass cannot promote alpha; unmet sample conditions remain inconclusive without extension.",
        },
    }


def verify_baselines_unchanged(contract: dict) -> dict:
    checks = {}
    for name, record in contract["baseline_artifact_hashes"].items():
        observed = sha256_file(resolve_workspace_path(record["path"]))
        checks[name] = {
            "expected_sha256": record["sha256"],
            "observed_sha256": observed,
            "unchanged": observed == record["sha256"],
        }
    if not all(item["unchanged"] for item in checks.values()):
        raise ValueError("a parent/baseline artifact changed during replication")
    return checks


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--combined-output", type=Path, default=DEFAULT_COMBINED)
    parser.add_argument("--verify-snapshot", action="store_true")
    args = parser.parse_args(argv)

    contract, parent, contract_sha256, start_ms, end_ms = load_and_verify_contract(
        args.contract
    )
    manifest, records = load_or_create_snapshot(
        args.contract, contract_sha256, args.snapshot
    )
    manifest_sha256 = sha256_file(args.snapshot)
    if args.verify_snapshot:
        print(json.dumps({
            "snapshot": str(args.snapshot),
            "sha256": manifest_sha256,
            "files": len(manifest["files"]),
            "verified": True,
        }, indent=2))
        return 0

    output_paths = {
        screen_id: resolve_workspace_path(screen["output"])
        for screen_id, screen in contract["screens"].items()
    }
    for path in [args.combined_output, *output_paths.values()]:
        if path.exists():
            raise FileExistsError(f"refusing to overwrite replication result: {path}")

    liquidation_module = load_module(
        "replication_parent_liquidation", LIQUIDATION_MODULE
    )
    deribit_module = load_module("replication_parent_deribit", DERIBIT_MODULE)
    liquidation = run_liquidation(
        parent, records, start_ms, end_ms, liquidation_module
    )
    deribit = run_deribit(parent, records, start_ms, end_ms, deribit_module)
    baseline_checks = verify_baselines_unchanged(contract)

    common = {
        "replication_contract": {
            "path": str(args.contract.relative_to(WORKSPACE_ROOT)),
            "sha256": contract_sha256,
            "created_at_utc": contract["created_at_utc"],
            "classification": contract["classification"],
        },
        "parent_contract": contract["parent_contract"],
        "snapshot_manifest": {
            "path": str(args.snapshot.relative_to(WORKSPACE_ROOT)),
            "sha256": manifest_sha256,
            "created_at_utc": manifest["created_at_utc"],
        },
        "replication_runner": manifest["replication_runner"],
        "claim_boundary": contract["claim_boundary"],
        "baseline_artifacts_unchanged": baseline_checks,
    }
    liquidation = {**common, **liquidation}
    deribit = {**common, **deribit}
    write_json_new(
        output_paths["liquidation_btc_eth_post_onset"], liquidation
    )
    write_json_new(output_paths["deribit_intraday_taker_gex"], deribit)

    combined = {
        "schema_version": 1,
        "campaign": contract["campaign"],
        "created_at_utc": now_utc(),
        **common,
        "window": contract["observation_window"],
        "screens": {
            "liquidation_btc_eth_post_onset": {
                "result": str(output_paths["liquidation_btc_eth_post_onset"].relative_to(WORKSPACE_ROOT)),
                "sha256": sha256_file(output_paths["liquidation_btc_eth_post_onset"]),
                "status": liquidation["truth"]["status"],
                "sample_all_met": liquidation["sample_conditions"]["all_met"],
            },
            "deribit_intraday_taker_gex": {
                "result": str(output_paths["deribit_intraday_taker_gex"].relative_to(WORKSPACE_ROOT)),
                "sha256": sha256_file(output_paths["deribit_intraday_taker_gex"]),
                "status": deribit["decision"]["verdict"],
                "sample_all_met": deribit["sample_conditions"]["all_met"],
            },
        },
        "decision": {
            "prior_liquidation_falsification_unchanged": True,
            "no_strategy_promoted": True,
            "no_orders_or_account_access": True,
            "no_private_data": True,
        },
    }
    write_json_new(args.combined_output, combined)
    print(json.dumps({
        "contract_sha256": contract_sha256,
        "snapshot_sha256": manifest_sha256,
        "liquidation": combined["screens"]["liquidation_btc_eth_post_onset"],
        "deribit": combined["screens"]["deribit_intraday_taker_gex"],
        "combined_output": str(args.combined_output),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
