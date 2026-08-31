#!/usr/bin/env python3
"""Run one feasibility-corrected chronological test on frozen public-data prefixes.

This runner exists because the prior 90-minute replication made both registered
sample stops unreachable. It preserves that trial and its outputs, rejects
structurally infeasible windows, validates stream coverage before metric code,
and emits complete alignment and timestamp provenance.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import replicate_accelerated_intraday as BASE


WORKSPACE_ROOT = BASE.WORKSPACE_ROOT
EXPLORATIONS = BASE.EXPLORATIONS
RUNNER_PATH = Path(__file__).resolve()
DEFAULT_CONTRACT = EXPLORATIONS / "accelerated-intraday-feasibility-corrected-contract.json"
DEFAULT_COMBINED = EXPLORATIONS / "accelerated-intraday-feasibility-corrected-results.json"
BLOCK_MS = 15 * 60 * 1000
CLASSIFICATION = (
    "outcome-blind retrospective feasibility-corrected chronological test; "
    "not prospective, not confirmation, and not a rescue test"
)


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def validate_window_feasibility(
    start_ms: int,
    end_ms: int,
    required_liquidation_blocks: int,
    required_deribit_transitions: int,
) -> dict[str, Any]:
    """Prove both parent sample stops are reachable before any outcome access."""
    if end_ms <= start_ms:
        raise ValueError("window end must be after start")
    if start_ms % BLOCK_MS or end_ms % BLOCK_MS:
        raise ValueError("window boundaries must be aligned to parent 15-minute blocks")
    duration_ms = end_ms - start_ms
    if duration_ms % BLOCK_MS:
        raise ValueError("window duration must contain whole parent blocks")

    window_blocks = duration_ms // BLOCK_MS
    complete_signal_blocks = max(0, window_blocks - 1)
    max_liquidation_blocks = complete_signal_blocks
    max_deribit_transitions = 2 * complete_signal_blocks
    report = {
        "block_ms": BLOCK_MS,
        "window_duration_ms": duration_ms,
        "window_minutes": duration_ms // 60_000,
        "window_blocks": window_blocks,
        "max_complete_liquidation_onset_blocks": max_liquidation_blocks,
        "required_complete_liquidation_onset_blocks": required_liquidation_blocks,
        "max_complete_deribit_asset_block_transitions": max_deribit_transitions,
        "required_complete_deribit_asset_block_transitions": required_deribit_transitions,
        "all_registered_stops_reachable": (
            max_liquidation_blocks >= required_liquidation_blocks
            and max_deribit_transitions >= required_deribit_transitions
        ),
    }
    if not report["all_registered_stops_reachable"]:
        raise ValueError(
            "registered sample stop is structurally unreachable: "
            f"liquidation blocks {max_liquidation_blocks}/{required_liquidation_blocks}; "
            f"Deribit transitions {max_deribit_transitions}/{required_deribit_transitions}"
        )
    return report


def _coverage(times: list[int], end_ms: int) -> dict[str, Any]:
    if not times:
        return {
            "observations": 0,
            "latest_ms": None,
            "latest_utc": None,
            "first_at_or_after_end_ms": None,
            "first_at_or_after_end_utc": None,
            "covers_end": False,
        }
    at_or_after = [timestamp for timestamp in times if timestamp >= end_ms]
    first_after = min(at_or_after) if at_or_after else None
    latest = max(times)
    return {
        "observations": len(times),
        "latest_ms": latest,
        "latest_utc": BASE.iso_ms(latest),
        "first_at_or_after_end_ms": first_after,
        "first_at_or_after_end_utc": BASE.iso_ms(first_after),
        "covers_end": first_after is not None,
    }


def validate_common_cutoffs(
    records: dict[str, list[dict]], end_ms: int
) -> dict[str, Any]:
    """Fail closed unless every sparse/event and BTC/ETH price stream covers end."""
    event_times: list[int] = []
    trade_times: list[int] = []
    index_times: dict[str, list[int]] = {"btc_usd": [], "eth_usd": []}

    for path, rows in records.items():
        name = Path(path).name
        if name.startswith("liquidation-capture-") or name == "liquidation-capture.jsonl":
            for row in rows:
                timestamp = BASE.force_order_event_ms(row)
                if _is_int(timestamp):
                    event_times.append(timestamp)
        elif name == "trades.jsonl":
            for row in rows:
                timestamp = (row.get("trade") or {}).get("timestamp")
                if _is_int(timestamp):
                    trade_times.append(timestamp)
        elif name == "index.jsonl":
            for row in rows:
                data = row.get("data") or {}
                index_name = data.get("index_name")
                timestamp = data.get("timestamp")
                if index_name in index_times and _is_int(timestamp):
                    index_times[index_name].append(timestamp)

    streams = {
        "liquidation_force_order": _coverage(event_times, end_ms),
        "deribit_option_trade": _coverage(trade_times, end_ms),
        "btc_usd_index": _coverage(index_times["btc_usd"], end_ms),
        "eth_usd_index": _coverage(index_times["eth_usd"], end_ms),
    }
    failed = [name for name, report in streams.items() if not report["covers_end"]]
    if failed:
        raise ValueError(
            "required stream does not establish coverage through exclusive end: "
            + ", ".join(failed)
        )
    common_cutoff_ms = min(report["latest_ms"] for report in streams.values())
    return {
        "required_end_ms": end_ms,
        "required_end_utc": BASE.iso_ms(end_ms),
        "common_cutoff_ms": common_cutoff_ms,
        "common_cutoff_utc": BASE.iso_ms(common_cutoff_ms),
        "all_required_streams_cover_end": True,
        "streams": streams,
    }


def alignment_provenance(outcomes: list[dict]) -> dict[str, Any]:
    """Serialize every incomplete onset and the exact absent alignment reason."""
    counts = {
        "initial": 0,
        "horizon_1m": 0,
        "horizon_5m": 0,
        "horizon_15m": 0,
    }
    incomplete: list[dict[str, Any]] = []
    for outcome in outcomes:
        if outcome["cohort_complete"]:
            continue
        alignments: dict[str, dict[str, Any]] = {}
        for label in ("0", "1m", "5m", "15m"):
            source = outcome[f"P_{label}"]
            alignments[f"P_{label}"] = {
                key: source[key]
                for key in (
                    "target_ms",
                    "observed_ts_ms",
                    "price",
                    "observed",
                    "missing",
                )
                if key in source
            }
            if "price" not in source:
                count_key = "initial" if label == "0" else f"horizon_{label}"
                counts[count_key] += 1
        incomplete.append(
            {
                "event_id": outcome["event_id"],
                "symbol": outcome["symbol"],
                "side": outcome["side"],
                "onset_ms": outcome["onset_ms"],
                "onset_utc": BASE.iso_ms(outcome["onset_ms"]),
                "alignments": alignments,
            }
        )
    return {
        "eligible_onsets": len(outcomes),
        "complete_onsets": sum(row["cohort_complete"] for row in outcomes),
        "incomplete_onsets_count": len(incomplete),
        "incomplete_alignment_counts": counts,
        "incomplete_onsets": incomplete,
    }


def instrument_clock_provenance(rows: list[dict]) -> dict[str, Any]:
    """Keep exchange creation time distinct from local capture availability."""
    creation_times = [
        (row.get("instrument") or {}).get("creation_timestamp") for row in rows
    ]
    creation_times = [value for value in creation_times if _is_int(value)]
    receipt_times = [row.get("received_ts_ms") for row in rows]
    receipt_times = [value for value in receipt_times if _is_int(value)]
    if not creation_times or not receipt_times:
        raise ValueError("instrument provenance requires creation and receipt clocks")
    return {
        "exchange_clock": {
            "field": "instrument.creation_timestamp",
            "first_ms": min(creation_times),
            "first_utc": BASE.iso_ms(min(creation_times)),
            "latest_ms": max(creation_times),
            "latest_utc": BASE.iso_ms(max(creation_times)),
        },
        "capture_availability_clock": {
            "field": "received_ts_ms",
            "first_ms": min(receipt_times),
            "first_utc": BASE.iso_ms(min(receipt_times)),
            "latest_ms": max(receipt_times),
            "latest_utc": BASE.iso_ms(max(receipt_times)),
        },
        "metadata_cutoff_policy": (
            "Exchange creation time is provenance. received_ts_ms is capture "
            "availability and determines no-lookahead metadata eligibility."
        ),
    }


def _resolve(record: dict) -> Path:
    return BASE.resolve_workspace_path(record["path"])


def _verify_hash_record(record: dict) -> None:
    observed = BASE.sha256_file(_resolve(record))
    if observed != record["sha256"]:
        raise ValueError(
            f"frozen file changed: {record['path']}: {observed} != {record['sha256']}"
        )


def load_and_verify_contract(path: Path) -> tuple[dict, dict, str, int, int, dict]:
    contract_bytes = path.read_bytes()
    contract_sha256 = BASE.sha256_bytes(contract_bytes)
    contract = json.loads(contract_bytes)
    if contract.get("status") != "frozen-before-corrected-outcome-analysis":
        raise ValueError("feasibility-corrected contract is not frozen")

    runner = contract["runner"]
    if _resolve(runner) != RUNNER_PATH or BASE.sha256_file(RUNNER_PATH) != runner["sha256"]:
        raise ValueError("corrected runner hash mismatch")
    for key in ("parent_contract", "source_snapshot", "review_audit"):
        _verify_hash_record(contract[key])
    for section in ("frozen_rule_files", "baseline_artifact_hashes"):
        for record in contract[section].values():
            _verify_hash_record(record)

    parent = json.loads(_resolve(contract["parent_contract"]).read_text(encoding="utf-8"))
    start_ms = BASE.parse_iso_ms(contract["observation_window"]["start_utc_inclusive"])
    end_ms = BASE.parse_iso_ms(contract["observation_window"]["end_utc_exclusive"])
    liquidation_blocks = parent["screens"]["liquidation_btc_eth_post_onset"]["sample_stop"][
        "minimum_nonempty_15m_onset_blocks"
    ]
    deribit_transitions = parent["screens"]["deribit_intraday_taker_gex"]["sample_stop"][
        "minimum_complete_asset_block_transitions"
    ]
    feasibility = validate_window_feasibility(
        start_ms, end_ms, liquidation_blocks, deribit_transitions
    )
    return contract, parent, contract_sha256, start_ms, end_ms, feasibility


def load_source_snapshot(contract: dict) -> tuple[dict, dict[str, list[dict]]]:
    snapshot_path = _resolve(contract["source_snapshot"])
    manifest = json.loads(snapshot_path.read_text(encoding="utf-8"))
    records = {
        item["path"]: BASE.read_frozen_snapshot(item) for item in manifest["files"]
    }
    if len(records) != contract["source_snapshot"]["file_count"]:
        raise ValueError("source snapshot file count mismatch")
    return manifest, records


def instrument_file_provenance(
    records: dict[str, list[dict]], manifest: dict
) -> list[dict[str, Any]]:
    manifest_by_path = {item["path"]: item for item in manifest["files"]}
    reports: list[dict[str, Any]] = []
    for path, rows in sorted(records.items()):
        if Path(path).name != "instruments.jsonl":
            continue
        source = manifest_by_path[path]
        reports.append(
            {
                "path": path,
                "prefix_bytes": source["prefix_bytes"],
                "sha256": source["sha256"],
                **instrument_clock_provenance(rows),
            }
        )
    if not reports:
        raise ValueError("source snapshot has no instrument metadata")
    return reports


def run_liquidation(
    parent: dict,
    records: dict[str, list[dict]],
    start_ms: int,
    end_ms: int,
    module: ModuleType,
) -> dict[str, Any]:
    event_rows = [
        row
        for path, rows in sorted(records.items())
        if Path(path).name.startswith("liquidation-capture-")
        for row in rows
    ]
    context_rows = BASE.filter_liquidation_context(event_rows, end_ms)
    index_rows, index_exclusions = BASE.filter_liquidation_index(
        BASE.records_by_filename(records, "index.jsonl"), start_ms, end_ms
    )
    cutoff_ms = end_ms - 1
    onset_stats = module.build_onsets(context_rows, cutoff_ms)
    eligible = BASE.filter_liquidation_onsets(
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
        "classification": CLASSIFICATION,
        "window": {
            "start_utc_inclusive": BASE.iso_ms(start_ms),
            "end_utc_exclusive": BASE.iso_ms(end_ms),
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
            "stop_utc": BASE.iso_ms(stop["stop_ms"]),
            "candidate_attempts": stop["tried"],
        },
        "sample_conditions": {
            "all_met": stop["found"],
            "conditions": conditions["conditions"],
            "by_symbol": conditions["by_symbol"],
            "by_side": conditions["by_side"],
            "nonempty_15m_blocks": conditions["nonempty_15m_blocks"],
        },
        "alignment_provenance": alignment_provenance(outcomes),
        "cohort": {
            "size": len(cohort),
            "event_ids": [row["event_id"] for row in cohort],
            "detail": cohort,
        },
        "gates": gates,
        "truth": truth,
        "decision_boundary": (
            "A pass cannot reverse the parent falsification or promote alpha. "
            "This protocol correction replaces only the structurally unreachable replication."
        ),
    }


def verify_baselines(contract: dict) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    for name, record in contract["baseline_artifact_hashes"].items():
        observed = BASE.sha256_file(_resolve(record))
        checks[name] = {
            "expected_sha256": record["sha256"],
            "observed_sha256": observed,
            "unchanged": observed == record["sha256"],
        }
    if not all(item["unchanged"] for item in checks.values()):
        raise ValueError("a frozen baseline artifact changed")
    return checks


def _preflight_outputs(paths: list[Path]) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(f"refusing to overwrite frozen artifact(s): {existing}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--combined-output", type=Path, default=DEFAULT_COMBINED)
    parser.add_argument("--verify-inputs", action="store_true")
    args = parser.parse_args(argv)

    (
        contract,
        parent,
        contract_sha256,
        start_ms,
        end_ms,
        feasibility,
    ) = load_and_verify_contract(args.contract)
    manifest, records = load_source_snapshot(contract)
    coverage = validate_common_cutoffs(records, end_ms)
    instrument_provenance = instrument_file_provenance(records, manifest)

    if args.verify_inputs:
        print(
            json.dumps(
                {
                    "contract": str(args.contract),
                    "contract_sha256": contract_sha256,
                    "source_snapshot_sha256": contract["source_snapshot"]["sha256"],
                    "files": len(records),
                    "feasibility": feasibility,
                    "coverage": coverage,
                    "verified": True,
                },
                indent=2,
            )
        )
        return 0

    output_paths = {
        screen: BASE.resolve_workspace_path(specification["output"])
        for screen, specification in contract["screens"].items()
    }
    _preflight_outputs([*output_paths.values(), args.combined_output])

    liquidation_module = BASE.load_module(
        "corrected_replication_parent_liquidation", BASE.LIQUIDATION_MODULE
    )
    deribit_module = BASE.load_module(
        "corrected_replication_parent_deribit", BASE.DERIBIT_MODULE
    )
    liquidation = run_liquidation(
        parent, records, start_ms, end_ms, liquidation_module
    )
    deribit = BASE.run_deribit(parent, records, start_ms, end_ms, deribit_module)
    deribit["classification"] = CLASSIFICATION
    baseline_checks = verify_baselines(contract)

    common = {
        "corrected_contract": {
            "path": str(args.contract.relative_to(WORKSPACE_ROOT)),
            "sha256": contract_sha256,
            "created_at_utc": contract["created_at_utc"],
            "classification": contract["classification"],
        },
        "parent_contract": contract["parent_contract"],
        "source_snapshot": contract["source_snapshot"],
        "review_audit": contract["review_audit"],
        "replication_runner": contract["runner"],
        "protocol_correction": contract["protocol_correction"],
        "claim_boundary": contract["claim_boundary"],
        "integrity_preconditions": {
            "window_feasibility": feasibility,
            "common_cutoff": coverage,
            "instrument_clock_provenance": instrument_provenance,
        },
        "baseline_artifacts_unchanged": baseline_checks,
    }
    liquidation = {**common, **liquidation}
    deribit = {**common, **deribit}

    BASE.write_json_new(
        output_paths["liquidation_btc_eth_post_onset"], liquidation
    )
    BASE.write_json_new(output_paths["deribit_intraday_taker_gex"], deribit)
    result_hashes = {
        "liquidation": BASE.sha256_file(
            output_paths["liquidation_btc_eth_post_onset"]
        ),
        "deribit": BASE.sha256_file(output_paths["deribit_intraday_taker_gex"]),
    }
    combined = {
        "schema_version": 1,
        "campaign": contract["campaign"],
        "created_at_utc": BASE.now_utc(),
        **common,
        "window": contract["observation_window"],
        "results": {
            "liquidation_btc_eth_post_onset": {
                "path": contract["screens"]["liquidation_btc_eth_post_onset"][
                    "output"
                ],
                "sha256": result_hashes["liquidation"],
                "status": liquidation["truth"]["status"],
                "sample_all_met": liquidation["sample_conditions"]["all_met"],
            },
            "deribit_intraday_taker_gex": {
                "path": contract["screens"]["deribit_intraday_taker_gex"]["output"],
                "sha256": result_hashes["deribit"],
                "status": deribit["decision"]["verdict"],
                "sample_all_met": deribit["sample_conditions"]["all_met"],
            },
        },
        "decision_boundary": contract["claim_boundary"],
    }
    BASE.write_json_new(args.combined_output, combined)
    print(
        json.dumps(
            {
                "contract_sha256": contract_sha256,
                "source_snapshot_sha256": contract["source_snapshot"]["sha256"],
                "liquidation": combined["results"][
                    "liquidation_btc_eth_post_onset"
                ],
                "deribit": combined["results"]["deribit_intraday_taker_gex"],
                "combined_output": str(args.combined_output),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
