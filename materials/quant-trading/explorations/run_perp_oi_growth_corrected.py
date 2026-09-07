#!/usr/bin/env python3
"""Frozen two-window test of cross-sectional open-interest growth on Binance USDT perpetuals,
with corrected boundary-row semantics (cycle 6).

Every rule except the handling of rows whose create_time falls outside the archive date is
imported from the frozen cycle-5 runner and pinned by hash; nothing is re-derived here. All
frozen hashes and upstream archive checks are verified before any candidate outcome is parsed.
Results use exclusive-create semantics. A pass is replicated; net-of-cost expected-return
evidence only, never deployable alpha.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import run_perp_cross_section_daily_reversal as parent
import run_perp_oi_growth_cross_section as oi5

WORKSPACE_ROOT = parent.WORKSPACE_ROOT
EXPLORATIONS = parent.EXPLORATIONS
CONTRACT_PATH = EXPLORATIONS / "pre_statement_perp_oi_growth_corrected.json"
INVENTORY_PATH = oi5.INVENTORY_PATH
KLINE_INVENTORY_PATH = oi5.KLINE_INVENTORY_PATH
AUDIT_PATH = EXPLORATIONS / "perp-oi-growth-corrected-prior-window-audit.json"
POLICY_PATH = parent.POLICY_PATH
PARENT_PATH = Path(parent.__file__).resolve()
OI5_PATH = Path(oi5.__file__).resolve()

DAY_MS = oi5.DAY_MS
HOUR_MS = oi5.HOUR_MS
ELIGIBILITY_DAYS = oi5.ELIGIBILITY_DAYS
END_OF_DAY_CUTOFF_MS = oi5.END_OF_DAY_CUTOFF_MS
MINIMUM_DATES = oi5.MINIMUM_DATES
MINIMUM_POSITIONS = oi5.MINIMUM_POSITIONS
METRICS_HEADER = oi5.METRICS_HEADER
REJECTION_KEYS = oi5.REJECTION_KEYS
ELIGIBILITY_KEYS = ("kline-history-incomplete", "below-liquidity-floor",
                    "missing-open-interest", "archive-without-end-of-day-row")
LOOP_ALPHA = 0.05 / 42
FROZEN_PARAMETERS = {**oi5.FROZEN_PARAMETERS, "loop_alpha": LOOP_ALPHA}
PIN_FIELDS = ("input_inventory_sha256", "kline_inventory_sha256", "prior_window_audit_sha256",
              "loop_policy_sha256", "parent_runner_sha256", "oi5_runner_sha256")

sha256_file = parent.sha256_file
parse_iso_ms = parent.parse_iso_ms
iso_ms = parent.iso_ms
write_json_new = parent.write_json_new
expected_dates = parent.expected_dates
block_bootstrap_p = parent.block_bootstrap_p
validate_primary_artifact = parent.validate_primary_artifact
_positive = parent._positive
_mean = parent._mean
_sd = parent._sd
_quarter_stability = parent._quarter_stability
_yearly_means = parent._yearly_means
_inside_root = parent._inside_root

_entry_symbol_date = oi5._entry_symbol_date
verify_metrics_archive = oi5.verify_metrics_archive
window_universe = oi5.window_universe
build_date = oi5.build_date
bootstrappable = oi5.bootstrappable


def check_contract_parameters(contract: dict) -> None:
    frozen = contract["frozen_rule"]["parameters"]
    if frozen != FROZEN_PARAMETERS:
        raise ValueError(
            f"contract parameters differ from runner constants: "
            f"contract {frozen} runner {FROZEN_PARAMETERS}")
    reused = {key: value for key, value in frozen.items() if key != "loop_alpha"}
    inherited = {key: value for key, value in oi5.FROZEN_PARAMETERS.items() if key != "loop_alpha"}
    if reused != inherited:
        raise ValueError("only loop_alpha may differ from the cycle-5 parameter set")


def evaluate_gates(rows: list[dict]) -> dict:
    """Cycle-5 gates with this cycle's alpha; G4 is strictly harder than at cycle 5."""
    gates = oi5.evaluate_gates(rows)
    g4 = gates["G4_dependence_robustness"]
    p_value = g4["one_sided_block_bootstrap_p_spread_net"]
    gates["G4_dependence_robustness"] = {**g4, "pass": p_value < LOOP_ALPHA, "alpha": LOOP_ALPHA}
    return gates


def verify_static_inputs(contract: dict) -> tuple[dict, dict, dict]:
    check_contract_parameters(contract)
    expected = {
        "input_inventory_sha256": contract["assets_and_data"]["input_inventory_sha256"],
        "kline_inventory_sha256": contract["assets_and_data"]["kline_inventory_sha256"],
        "prior_window_audit_sha256": contract["novelty_and_consumption"]["prior_window_audit_sha256"],
        "loop_policy_sha256": contract["loop_multiplicity"]["policy_sha256"],
        "parent_runner_sha256": contract["governance"]["parent_runner_sha256"],
        "oi5_runner_sha256": contract["governance"]["oi5_runner_sha256"],
    }
    observed = {
        "input_inventory_sha256": sha256_file(INVENTORY_PATH),
        "kline_inventory_sha256": sha256_file(KLINE_INVENTORY_PATH),
        "prior_window_audit_sha256": sha256_file(AUDIT_PATH),
        "loop_policy_sha256": sha256_file(POLICY_PATH),
        "parent_runner_sha256": sha256_file(PARENT_PATH),
        "oi5_runner_sha256": sha256_file(OI5_PATH),
    }
    for field, value in expected.items():
        if observed[field] != value:
            raise ValueError(f"{field} mismatch: expected {value}, observed {observed[field]}")

    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
    if inventory.get("candidate_outcomes_accessed") is not False:
        raise ValueError("input inventory does not certify pre-outcome construction")
    if inventory.get("kline_inventory_sha256") != observed["kline_inventory_sha256"]:
        raise ValueError("input inventory does not bind the frozen kline inventory")
    windows = inventory.get("windows", {})
    if set(windows) != {"primary", "replication"}:
        raise ValueError("input inventory must carry exactly the primary and replication universes")
    for label, window in windows.items():
        universe = window.get("universe", [])
        if len(universe) != oi5.UNIVERSE_SIZE or len(set(universe)) != oi5.UNIVERSE_SIZE:
            raise ValueError(f"{label} universe must hold {oi5.UNIVERSE_SIZE} distinct symbols")
    entries = inventory.get("files", [])
    if inventory.get("archive_count") != len(entries) or not entries:
        raise ValueError("input inventory archive count mismatch")
    paths = [entry["path"] for entry in entries]
    if len(paths) != len(set(paths)):
        raise ValueError("input inventory contains duplicate paths")
    allowed = set()
    for window in windows.values():
        allowed.update(window["universe"])
    for entry in entries:
        symbol, _ = _entry_symbol_date(entry)
        if symbol not in allowed:
            raise ValueError(f"inventory archive outside both universes: {entry['path']}")
        verify_metrics_archive(entry)

    kline_inventory = json.loads(KLINE_INVENTORY_PATH.read_text(encoding="utf-8"))
    if kline_inventory.get("candidate_outcomes_accessed") is not False:
        raise ValueError("kline inventory does not certify pre-outcome construction")
    for entry in kline_inventory["files"]:
        parent.verify_archive(entry)
    return inventory, kline_inventory, {
        **observed, "metrics_archives_verified": len(entries),
        "kline_archives_verified": len(kline_inventory["files"])}


def _parse_stamp(stamp: str) -> tuple[str, int]:
    """(date part, milliseconds into that date). Malformed or out-of-range stamps abort.

    Strictly stricter than the cycle-5 parser: every digit position must be an ASCII digit, so
    space- or sign-padded fields are refused instead of silently coerced. This parser never
    selects a value the cycle-5 parser would not have selected; wherever the two disagree this
    one aborts and never substitutes. An outcome-blind scan of every create_time string in all
    51,224 frozen archives found zero padded or non-canonical stamps, so the check is inert here.
    """
    if len(stamp) != 19 or stamp[10] != " " or stamp[13] != ":" or stamp[16] != ":":
        raise ValueError(f"malformed create_time: {stamp!r}")
    if stamp[4] != "-" or stamp[7] != "-":
        raise ValueError(f"malformed create_time date: {stamp!r}")
    digits = stamp[:4] + stamp[5:7] + stamp[8:10] + stamp[11:13] + stamp[14:16] + stamp[17:19]
    if not digits.isdigit() or not digits.isascii():
        raise ValueError(f"non-digit create_time field: {stamp!r}")
    date, hour, minute, second = stamp[:10], int(stamp[11:13]), int(stamp[14:16]), int(stamp[17:19])
    datetime.strptime(date, "%Y-%m-%d")
    if not (0 <= hour < 24 and 0 <= minute < 60 and 0 <= second < 60):
        raise ValueError(f"create_time fields out of range: {stamp!r}")
    return date, (hour * 3600 + minute * 60 + second) * 1000


def _parse_metrics_archive(entry: dict, census: dict[str, int],
                           root: Path = WORKSPACE_ROOT) -> float | None:
    """End-of-day open interest: the last in-date row at or after the cutoff.

    A row whose create_time date differs from the archive date is dropped and counted; it can
    never be selected as end of day. Nothing is imputed. Malformed stamps, symbol mismatches,
    column-count drift, and header drift still abort the invocation.
    """
    symbol, date = _entry_symbol_date(entry)
    last_ms = -1
    last_value: float | None = None
    off_date = 0
    with zipfile.ZipFile(root / entry["path"]) as archive:
        with archive.open(entry["zip_member"]) as raw:
            reader = csv.reader(io.TextIOWrapper(raw, encoding="utf-8"))
            header = next(reader)
            if tuple(header) != METRICS_HEADER:
                raise ValueError(f"metrics header drift: {entry['path']}")
            for line_number, row in enumerate(reader, 2):
                if len(row) != len(METRICS_HEADER):
                    raise ValueError(f"schema column count at {entry['path']}:{line_number}")
                if row[1] != symbol:
                    raise ValueError(f"symbol column mismatch at {entry['path']}:{line_number}")
                try:
                    row_date, time_ms = _parse_stamp(row[0])
                    value = float(row[2])
                except ValueError as exc:
                    raise ValueError(f"schema value at {entry['path']}:{line_number}: {exc}") from exc
                if row_date != date:
                    off_date += 1
                    census["off-date-boundary-row"] += 1
                    offset = (datetime.strptime(row_date, "%Y-%m-%d").date()
                              - datetime.strptime(date, "%Y-%m-%d").date()).days
                    census[f"off-date-offset-days:{offset:+d}"] = (
                        census.get(f"off-date-offset-days:{offset:+d}", 0) + 1)
                    if time_ms >= END_OF_DAY_CUTOFF_MS:
                        census["off-date-row-at-or-after-cutoff"] += 1
                    continue
                if time_ms >= END_OF_DAY_CUTOFF_MS and time_ms >= last_ms:
                    last_ms = time_ms
                    last_value = value
    if off_date:
        census["archives-with-off-date-row"] += 1
    if last_ms < 0:
        return None
    return _positive(last_value)


def load_end_of_day_open_interest(inventory: dict, symbols: list[str], lo_ms: int, hi_ms: int,
                                  counters: dict[str, int],
                                  census: dict[str, int]) -> dict[tuple[str, int], float]:
    wanted = set(symbols)
    seen: set[tuple[str, str]] = set()
    values: dict[tuple[str, int], float] = {}
    for entry in inventory["files"]:
        symbol, date = _entry_symbol_date(entry)
        if symbol not in wanted:
            continue
        day_ms = parse_iso_ms(date + "T00:00:00Z")
        if not lo_ms <= day_ms < hi_ms:
            continue
        if (symbol, date) in seen:
            raise ValueError(f"duplicate inventory date for {symbol} {date}")
        seen.add((symbol, date))
        census["archives-read"] += 1
        value = _parse_metrics_archive(entry, census)
        if value is None:
            counters["archive-without-end-of-day-row"] += 1
            continue
        values[(symbol, day_ms)] = value
    return values


def run_window(contract: dict, label: str, inventory: dict, kline_inventory: dict) -> dict:
    window = contract["windows"][label]
    start_ms = parse_iso_ms(window["start_utc_inclusive"])
    end_ms = parse_iso_ms(window["end_utc_exclusive"])
    dates = expected_dates(start_ms, end_ms)
    registered = inventory["windows"][label]
    if (registered["start_utc_inclusive"], registered["end_utc_exclusive"]) != (
            window["start_utc_inclusive"], window["end_utc_exclusive"]):
        raise ValueError(f"{label} window differs between contract and input inventory")
    universe, lookback = window_universe(kline_inventory, start_ms)
    if universe != registered["universe"]:
        raise ValueError(f"{label} universe recomputed from frozen klines differs from the inventory")
    counters = {key: 0 for key in ELIGIBILITY_KEYS}
    census = {key: 0 for key in ("archives-read", "archives-with-off-date-row",
                                 "off-date-boundary-row", "off-date-row-at-or-after-cutoff")}
    rejected = {key: 0 for key in REJECTION_KEYS}
    records = parent.load_day_records(kline_inventory, start_ms - ELIGIBILITY_DAYS * DAY_MS,
                                      end_ms + 2 * HOUR_MS)
    open_interest = load_end_of_day_open_interest(inventory, universe, start_ms - 2 * DAY_MS,
                                                  end_ms, counters, census)
    rows: list[dict] = []
    for day_ms in dates:
        row = build_date(records, open_interest, universe, day_ms, counters, rejected)
        if row is not None:
            rows.append(row)
    positions = sum(row["positions"] for row in rows)
    sample = {
        "expected_dates": len(dates), "included_dates": len(rows), "positions": positions,
        "positions_by_side": {"long": sum(row["long_positions"] for row in rows),
                              "short": sum(row["short_positions"] for row in rows)},
        "included_quarters": len({row["quarter"] for row in rows}),
        "mean_eligible": _mean([row["eligible"] for row in rows]) if rows else 0.0,
        "mean_leg_size": _mean([row["leg_size"] for row in rows]) if rows else 0.0,
        "distinct_symbols_held": len({symbol for row in rows for symbol in row["symbols"]}),
        "end_of_day_open_interest_observations": len(open_interest),
    }
    sample_met = len(rows) >= MINIMUM_DATES and positions >= MINIMUM_POSITIONS
    if not sample_met:
        gates: dict = {}
        status = "inconclusive"
    else:
        gates = evaluate_gates(rows)
        status = ("pass-on-this-sample" if all(gate["pass"] for gate in gates.values())
                  else "falsified")
    spread_values = [row["spread_net"] for row in rows]
    lag_values = [row["lag_spread_net"] for row in rows if row["lag_spread_net"] is not None]
    pooled = {
        "mean_spread_net": _mean(spread_values) if rows else None,
        "mean_spread_gross": _mean([row["spread_gross"] for row in rows]) if rows else None,
        "mean_long_net": _mean([row["long_net"] for row in rows]) if rows else None,
        "mean_short_net": _mean([row["short_net"] for row in rows]) if rows else None,
        "mean_long_excess": _mean([row["long_excess"] for row in rows]) if rows else None,
        "mean_short_excess": _mean([row["short_excess"] for row in rows]) if rows else None,
        "mean_market_gross": _mean([row["market_gross"] for row in rows]) if rows else None,
        "spread_net_sd": _sd(spread_values),
    }
    diagnostics = {
        "quarterly_means": _quarter_stability(rows, "spread_net")["quarterly_means"] if rows else {},
        "yearly_mean_spread_net": _yearly_means(rows, "spread_net"),
        "mean_leg_gross": {"long": _mean([row["long_gross"] for row in rows]) if rows else None,
                           "short": _mean([row["short_gross"] for row in rows]) if rows else None},
        "gross_spread_block_bootstrap_p": block_bootstrap_p(
            {row["date"]: row["spread_gross"] for row in rows},
            replicates=oi5.BOOTSTRAP_REPLICATES, seed=oi5.BOOTSTRAP_SEED,
            block_dates=oi5.BOOTSTRAP_BLOCK_DATES) if bootstrappable(rows) else None,
        "one_hour_lagged_net_spread": {
            "dates": len(lag_values), "mean": _mean(lag_values) if lag_values else None,
            "one_sided_block_bootstrap_p": block_bootstrap_p(
                {row["date"]: row["lag_spread_net"] for row in rows if row["lag_spread_net"] is not None},
                replicates=oi5.BOOTSTRAP_REPLICATES, seed=oi5.BOOTSTRAP_SEED,
                block_dates=oi5.BOOTSTRAP_BLOCK_DATES) if bootstrappable(lag_values) else None},
        "eligibility_exclusions": counters, "rejected": rejected,
        "archive_census": census,
        "funding_timestamp_crossings": {
            "per_position_upper_bound": oi5.FUNDING_CROSSINGS_PER_POSITION,
            "total_upper_bound": oi5.FUNDING_CROSSINGS_PER_POSITION * positions,
            "subtracted": False},
        "universe": universe, "universe_lookback": lookback,
    }
    return {
        "role": label, "window": [window["start_utc_inclusive"], window["end_utc_exclusive"]],
        "status": status, "sample_minimum_met": sample_met, "sample": sample, "gate_results": gates,
        "pooled": pooled, "diagnostics": diagnostics, "claim_boundary": contract["claim_boundary"],
        "dates": [{key: value for key, value in row.items() if key != "symbols"} for row in rows],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--window", choices=("primary", "replication"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--primary-artifact", type=Path)
    parser.add_argument("--primary-artifact-sha256")
    parser.add_argument("--verify-inputs-only", action="store_true")
    args = parser.parse_args(argv)

    observed_contract = sha256_file(CONTRACT_PATH)
    if observed_contract != args.contract_sha256:
        raise SystemExit(f"contract hash mismatch: expected {args.contract_sha256}, observed {observed_contract}")
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"refusing to overwrite existing result: {output}")
    _inside_root(output, WORKSPACE_ROOT)
    inventory, kline_inventory, verified = verify_static_inputs(contract)
    runner_sha = sha256_file(Path(__file__).resolve())
    pins = {"runner_sha256": runner_sha, **{field: verified[field] for field in PIN_FIELDS}}
    if args.verify_inputs_only:
        print(json.dumps({**pins, "contract_sha256": observed_contract,
                          "metrics_archives_verified": verified["metrics_archives_verified"],
                          "kline_archives_verified": verified["kline_archives_verified"]}, indent=2, sort_keys=True))
        return 0

    primary_sha = None
    if args.window == "replication":
        if args.primary_artifact is None or not args.primary_artifact_sha256:
            raise SystemExit("replication requires --primary-artifact and --primary-artifact-sha256")
        primary_sha = sha256_file(args.primary_artifact)
        if primary_sha != args.primary_artifact_sha256:
            raise SystemExit("primary artifact hash mismatch")
        primary = json.loads(args.primary_artifact.read_text(encoding="utf-8"))
        validate_primary_artifact(primary, {"contract_sha256": observed_contract, "runner_sha256": runner_sha})
        primary_pins = {field: primary.get("static_inputs", {}).get(field) for field in PIN_FIELDS}
        if primary_pins != {field: verified[field] for field in PIN_FIELDS}:
            raise SystemExit(f"primary artifact static-input pins differ from this invocation: "
                             f"primary {primary_pins} now {pins}")
    result = run_window(contract, args.window, inventory, kline_inventory)
    payload = {
        "schema_version": 1, "contract_id": contract["id"], "contract_sha256": observed_contract,
        "runner_sha256": runner_sha, "static_inputs": verified,
        "generated_at_utc": iso_ms(int(datetime.now(timezone.utc).timestamp() * 1000)),
        "role": args.window, "result": result,
    }
    if args.window == "replication":
        payload["primary_artifact"] = str(args.primary_artifact.resolve().relative_to(WORKSPACE_ROOT))
        payload["primary_artifact_sha256"] = primary_sha
    write_json_new(args.output, payload)
    print(json.dumps({f"{args.window}_status": result["status"], "included_dates": result["sample"]["included_dates"],
                      "positions": result["sample"]["positions"], "output": str(args.output)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
