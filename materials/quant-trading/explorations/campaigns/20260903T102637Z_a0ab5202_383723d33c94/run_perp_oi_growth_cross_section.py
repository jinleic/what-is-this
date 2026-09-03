#!/usr/bin/env python3
"""Frozen two-window test of cross-sectional open-interest growth on Binance USDT perpetuals.

Signal: prior-day log change in end-of-day open interest (contracts) from the
Data Vision USD-M daily 5-minute metrics archives. Long the top quintile, short
the bottom quintile among point-in-time eligible symbols of a pre-registered
72-symbol universe; enter at the UTC-day open, exit at the next UTC-day open,
net of two 0.05% taker fees per position. Kline parsing, the exit chain, fee
algebra, and the block bootstrap are reused from the hash-pinned cycle-4 runner.
All frozen hashes and every archive are verified before any outcome is read.
Results use exclusive-create semantics.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import run_perp_cross_section_daily_reversal as parent

WORKSPACE_ROOT = parent.WORKSPACE_ROOT
EXPLORATIONS = parent.EXPLORATIONS
CONTRACT_PATH = EXPLORATIONS / "pre_statement_perp_oi_growth.json"
INVENTORY_PATH = EXPLORATIONS / "perp-oi-growth-inputs.json"
KLINE_INVENTORY_PATH = parent.INVENTORY_PATH
AUDIT_PATH = EXPLORATIONS / "perp-oi-growth-prior-window-audit.json"
POLICY_PATH = parent.POLICY_PATH
PARENT_PATH = Path(parent.__file__).resolve()

HOUR_MS = parent.HOUR_MS
DAY_MS = parent.DAY_MS
TAKER_FEE = parent.TAKER_FEE
ELIGIBILITY_DAYS = 20
LIQUIDITY_FLOOR_USD = 5_000_000.0
UNIVERSE_SIZE = 72
MIN_COMPLETE_LOOKBACK_DAYS = 300
MINIMUM_CROSS_SECTION = 30
QUINTILE_DENOMINATOR = 5
END_OF_DAY_CUTOFF_MS = 23 * HOUR_MS + 45 * 60_000
MAGNITUDE_THRESHOLD = 0.0005
QUARTER_STABILITY_NUMERATOR = 3
QUARTER_STABILITY_DENOMINATOR = 5
MINIMUM_DATES = 300
MINIMUM_POSITIONS = 2_000
FUNDING_CROSSINGS_PER_POSITION = 3
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_SEED = 20_260_902
BOOTSTRAP_BLOCK_DATES = 7
LOOP_ALPHA = 0.05 / 30
METRICS_HEADER = (
    "create_time", "symbol", "sum_open_interest", "sum_open_interest_value",
    "count_toptrader_long_short_ratio", "sum_toptrader_long_short_ratio",
    "count_long_short_ratio", "sum_taker_long_short_vol_ratio",
)
FROZEN_PARAMETERS = {
    "taker_fee_per_side": TAKER_FEE,
    "eligibility_days": ELIGIBILITY_DAYS,
    "liquidity_floor_usd": LIQUIDITY_FLOOR_USD,
    "universe_size": UNIVERSE_SIZE,
    "minimum_complete_lookback_days": MIN_COMPLETE_LOOKBACK_DAYS,
    "minimum_cross_section": MINIMUM_CROSS_SECTION,
    "quintile_denominator": QUINTILE_DENOMINATOR,
    "end_of_day_cutoff_ms": END_OF_DAY_CUTOFF_MS,
    "magnitude_threshold": MAGNITUDE_THRESHOLD,
    "quarter_stability_ratio": [QUARTER_STABILITY_NUMERATOR, QUARTER_STABILITY_DENOMINATOR],
    "minimum_dates": MINIMUM_DATES,
    "minimum_positions": MINIMUM_POSITIONS,
    "funding_crossings_per_position": FUNDING_CROSSINGS_PER_POSITION,
    "bootstrap_replicates": BOOTSTRAP_REPLICATES,
    "bootstrap_seed": BOOTSTRAP_SEED,
    "bootstrap_block_dates": BOOTSTRAP_BLOCK_DATES,
    "loop_alpha": LOOP_ALPHA,
}
REJECTION_KEYS = ("insufficient-cross-section", "missing-entry-bar", "unvaluable-exit",
                  "early-exit-last-bar", "lag-missing-entry-bar", "lag-unvaluable-exit",
                  "lag-early-exit-last-bar", "lag-exit-at-next-midnight")
PIN_FIELDS = ("input_inventory_sha256", "kline_inventory_sha256", "prior_window_audit_sha256",
              "loop_policy_sha256", "parent_runner_sha256")

sha256_file = parent.sha256_file
parse_iso_ms = parent.parse_iso_ms
iso_ms = parent.iso_ms
date_label = parent.date_label
write_json_new = parent.write_json_new
expected_dates = parent.expected_dates
quarter_label = parent.quarter_label
is_complete = parent.is_complete
exit_price = parent.exit_price
net_return = parent.net_return
block_bootstrap_p = parent.block_bootstrap_p
validate_primary_artifact = parent.validate_primary_artifact
_positive = parent._positive
_mean = parent._mean
_sd = parent._sd
_quarter_stability = parent._quarter_stability
_yearly_means = parent._yearly_means
_inside_root = parent._inside_root


def check_contract_parameters(contract: dict) -> None:
    frozen = contract["frozen_rule"]["parameters"]
    if frozen != FROZEN_PARAMETERS:
        raise ValueError(
            f"contract parameters differ from runner constants: "
            f"contract {frozen} runner {FROZEN_PARAMETERS}")


def universe_lookback(start_ms: int) -> tuple[int, int]:
    """[1 January of the preceding calendar year, window start)."""
    start = datetime.fromtimestamp(start_ms / 1000, timezone.utc)
    if (start.month, start.day, start.hour, start.minute, start.second) != (1, 1, 0, 0, 0):
        raise ValueError("window start must be 1 January 00:00:00 UTC")
    lookback_start = int(start.replace(year=start.year - 1).timestamp() * 1000)
    return lookback_start, start_ms


def window_universe(kline_inventory: dict, start_ms: int) -> tuple[list[str], dict[str, dict]]:
    """Top UNIVERSE_SIZE symbols by preceding-calendar-year quote volume over complete days."""
    lookback_start, lookback_end = universe_lookback(start_ms)
    records = parent.load_day_records(kline_inventory, lookback_start, lookback_end)
    stats: dict[str, dict] = {}
    for (symbol, _day_ms), record in records.items():
        if is_complete(record):
            entry = stats.setdefault(symbol, {"complete_days": 0, "quote_volume_usd": 0.0})
            entry["complete_days"] += 1
            entry["quote_volume_usd"] += record["quote_volume"]
    ranked = sorted(
        (symbol for symbol, entry in stats.items()
         if entry["complete_days"] >= MIN_COMPLETE_LOOKBACK_DAYS),
        key=lambda symbol: (-stats[symbol]["quote_volume_usd"], symbol))
    chosen = ranked[:UNIVERSE_SIZE]
    return chosen, {symbol: stats[symbol] for symbol in chosen}


def _entry_symbol_date(entry: dict) -> tuple[str, str]:
    name = Path(entry["path"]).name
    if not name.endswith(".zip") or "-metrics-" not in name:
        raise ValueError(f"unexpected metrics archive name: {entry['path']}")
    symbol, suffix = name.split("-metrics-", 1)
    date = suffix[:-4]
    if symbol != entry["symbol"] or date != entry["date"]:
        raise ValueError(f"inventory symbol/date disagree with archive name: {entry['path']}")
    return symbol, date


def verify_metrics_archive(entry: dict, root: Path = WORKSPACE_ROOT) -> str:
    _entry_symbol_date(entry)
    archive_path = _inside_root(root / entry["path"], root)
    sidecar_path = _inside_root(root / entry["checksum_sidecar"], root)
    if not archive_path.is_file() or not sidecar_path.is_file():
        raise ValueError(f"archive or sidecar absent: {entry['path']}")
    if archive_path.stat().st_size != entry["size_bytes"]:
        raise ValueError(f"archive byte size drift: {entry['path']}")
    if sha256_file(sidecar_path) != entry["checksum_sidecar_sha256"]:
        raise ValueError(f"sidecar file hash drift: {entry['checksum_sidecar']}")
    fields = sidecar_path.read_text(encoding="utf-8").split()
    if len(fields) < 2:
        raise ValueError(f"malformed upstream sidecar: {entry['checksum_sidecar']}")
    if fields[0].lower() != entry["sha256"]:
        raise ValueError(f"upstream sidecar hash mismatch: {entry['path']}")
    if Path(fields[-1]).name != archive_path.name:
        raise ValueError(f"upstream sidecar filename mismatch: {entry['path']}")
    if sha256_file(archive_path) != entry["sha256"]:
        raise ValueError(f"archive hash drift: {entry['path']}")
    with zipfile.ZipFile(archive_path) as archive:
        members = archive.namelist()
        if members != [entry["zip_member"]]:
            raise ValueError(f"zip member mismatch: {entry['path']}")
        corrupt = archive.testzip()
        if corrupt is not None:
            raise ValueError(f"zip CRC failure: {entry['path']}:{corrupt}")
        with archive.open(members[0]) as raw:
            first = next(csv.reader(io.TextIOWrapper(raw, encoding="utf-8")))
    if tuple(first) != METRICS_HEADER:
        raise ValueError(f"metrics header drift: {entry['path']}")
    return members[0]


def verify_static_inputs(contract: dict) -> tuple[dict, dict, dict]:
    check_contract_parameters(contract)
    expected = {
        "input_inventory_sha256": contract["assets_and_data"]["input_inventory_sha256"],
        "kline_inventory_sha256": contract["assets_and_data"]["kline_inventory_sha256"],
        "prior_window_audit_sha256": contract["novelty_and_consumption"]["prior_window_audit_sha256"],
        "loop_policy_sha256": contract["loop_multiplicity"]["policy_sha256"],
        "parent_runner_sha256": contract["governance"]["parent_runner_sha256"],
    }
    observed = {
        "input_inventory_sha256": sha256_file(INVENTORY_PATH),
        "kline_inventory_sha256": sha256_file(KLINE_INVENTORY_PATH),
        "prior_window_audit_sha256": sha256_file(AUDIT_PATH),
        "loop_policy_sha256": sha256_file(POLICY_PATH),
        "parent_runner_sha256": sha256_file(PARENT_PATH),
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
        if len(universe) != UNIVERSE_SIZE or len(set(universe)) != UNIVERSE_SIZE:
            raise ValueError(f"{label} universe must hold {UNIVERSE_SIZE} distinct symbols")
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


def _parse_time_of_day(stamp: str, date: str) -> int:
    if len(stamp) != 19 or stamp[:10] != date or stamp[10] != " " or stamp[13] != ":" or stamp[16] != ":":
        raise ValueError(f"create_time outside archive date or malformed: {stamp!r}")
    hour, minute, second = int(stamp[11:13]), int(stamp[14:16]), int(stamp[17:19])
    if not (0 <= hour < 24 and 0 <= minute < 60 and 0 <= second < 60):
        raise ValueError(f"create_time fields out of range: {stamp!r}")
    return (hour * 3600 + minute * 60 + second) * 1000


def _parse_metrics_archive(entry: dict, root: Path = WORKSPACE_ROOT) -> float | None:
    """End-of-day open interest: the last row at or after the cutoff; None if absent or nonpositive."""
    symbol, date = _entry_symbol_date(entry)
    last_ms = -1
    last_value: float | None = None
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
                    time_ms = _parse_time_of_day(row[0], date)
                    value = float(row[2])
                except ValueError as exc:
                    raise ValueError(f"schema value at {entry['path']}:{line_number}: {exc}") from exc
                if time_ms >= END_OF_DAY_CUTOFF_MS and time_ms >= last_ms:
                    last_ms = time_ms
                    last_value = value
    if last_ms < 0:
        return None
    return _positive(last_value)


def load_end_of_day_open_interest(inventory: dict, symbols: list[str], lo_ms: int, hi_ms: int,
                                  counters: dict[str, int]) -> dict[tuple[str, int], float]:
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
        value = _parse_metrics_archive(entry)
        if value is None:
            counters["archive-without-end-of-day-row"] += 1
            continue
        values[(symbol, day_ms)] = value
    return values


def growth_signal(open_interest: dict, symbol: str, day_ms: int) -> float | None:
    """Prior-day log change in end-of-day open interest, known at the day-d open."""
    previous = open_interest.get((symbol, day_ms - DAY_MS))
    before = open_interest.get((symbol, day_ms - 2 * DAY_MS))
    if previous is None or before is None:
        return None
    return math.log(previous / before)


def eligible_symbols(records: dict, open_interest: dict, symbols: list[str], day_ms: int,
                     counters: dict[str, int]) -> dict[str, float]:
    """Point-in-time eligible symbols at 00:00 UTC of day_ms with their signal."""
    eligible: dict[str, float] = {}
    for symbol in symbols:
        total_quote = 0.0
        complete = True
        for back in range(1, ELIGIBILITY_DAYS + 1):
            record = records.get((symbol, day_ms - back * DAY_MS))
            if not is_complete(record):
                complete = False
                break
            total_quote += record["quote_volume"]
        if not complete:
            counters["kline-history-incomplete"] += 1
            continue
        if total_quote / ELIGIBILITY_DAYS < LIQUIDITY_FLOOR_USD:
            counters["below-liquidity-floor"] += 1
            continue
        signal = growth_signal(open_interest, symbol, day_ms)
        if signal is None:
            counters["missing-open-interest"] += 1
            continue
        eligible[symbol] = signal
    return eligible


def _leg_row(records: dict, symbol: str, day_ms: int, side: str, signal: float,
             rejected: dict[str, int]) -> dict | None:
    today = records.get((symbol, day_ms))
    entry = _positive(today["open0"]) if today is not None else None
    if entry is None:
        rejected["missing-entry-bar"] += 1
        return None
    exit_value, fallback = exit_price(records, symbol, day_ms)
    if exit_value is None:
        rejected["unvaluable-exit"] += 1
        return None
    if fallback is not None:
        rejected[fallback] += 1
    gross = exit_value / entry - 1.0
    row = {"symbol": symbol, "side": side, "signal": signal, "gross": gross,
           "net": net_return(gross, side), "exit_fallback": fallback}
    lagged_entry = _positive(today["open1"])
    if lagged_entry is None:
        rejected["lag-missing-entry-bar"] += 1
        return row
    lagged_exit, lagged_fallback = exit_price(records, symbol, day_ms, lagged=True)
    if lagged_exit is None:
        rejected["lag-unvaluable-exit"] += 1
        return row
    if lagged_fallback is not None:
        rejected[lagged_fallback] += 1
    row["lag_net"] = net_return(lagged_exit / lagged_entry - 1.0, side)
    return row


def build_date(records: dict, open_interest: dict, symbols: list[str], day_ms: int,
               counters: dict[str, int], rejected: dict[str, int]) -> dict | None:
    signals = eligible_symbols(records, open_interest, symbols, day_ms, counters)
    if len(signals) < MINIMUM_CROSS_SECTION:
        rejected["insufficient-cross-section"] += 1
        return None
    ranked = sorted(signals, key=lambda symbol: (signals[symbol], symbol))
    leg = len(ranked) // QUINTILE_DENOMINATOR
    positions: list[dict] = []
    for symbol in ranked[-leg:]:
        row = _leg_row(records, symbol, day_ms, "long", signals[symbol], rejected)
        if row is not None:
            positions.append(row)
    for symbol in ranked[:leg]:
        row = _leg_row(records, symbol, day_ms, "short", signals[symbol], rejected)
        if row is not None:
            positions.append(row)
    longs = [row for row in positions if row["side"] == "long"]
    shorts = [row for row in positions if row["side"] == "short"]
    if not longs or not shorts:
        rejected["empty-leg"] = rejected.get("empty-leg", 0) + 1
        return None
    market: list[float] = []
    for symbol in ranked:
        today = records.get((symbol, day_ms))
        entry = _positive(today["open0"]) if today is not None else None
        if entry is None:
            continue
        exit_value, _ = exit_price(records, symbol, day_ms)
        if exit_value is not None:
            market.append(exit_value / entry - 1.0)
    market_gross = _mean(market)
    long_net = _mean([row["net"] for row in longs])
    short_net = _mean([row["net"] for row in shorts])
    long_gross = _mean([row["gross"] for row in longs])
    short_gross = _mean([-row["gross"] for row in shorts])
    lag_longs = [row["lag_net"] for row in longs if "lag_net" in row]
    lag_shorts = [row["lag_net"] for row in shorts if "lag_net" in row]
    lag_spread = (_mean(lag_longs) + _mean(lag_shorts)) / 2.0 if lag_longs and lag_shorts else None
    date = date_label(day_ms)
    return {
        "date": date, "quarter": quarter_label(date), "eligible": len(ranked), "leg_size": leg,
        "positions": len(positions), "long_positions": len(longs), "short_positions": len(shorts),
        "market_gross": market_gross, "long_net": long_net, "short_net": short_net,
        "long_gross": long_gross, "short_gross": short_gross,
        "spread_net": (long_net + short_net) / 2.0, "spread_gross": (long_gross + short_gross) / 2.0,
        "long_excess": long_net - market_gross, "short_excess": short_net + market_gross,
        "lag_spread_net": lag_spread, "symbols": [row["symbol"] for row in positions],
    }


def bootstrappable(values: list) -> bool:
    return len(values) >= BOOTSTRAP_BLOCK_DATES * 3


def evaluate_gates(rows: list[dict]) -> dict:
    spread = {row["date"]: row["spread_net"] for row in rows}
    stability = _quarter_stability(rows, "spread_net")
    p_value = block_bootstrap_p(spread, replicates=BOOTSTRAP_REPLICATES, seed=BOOTSTRAP_SEED,
                                block_dates=BOOTSTRAP_BLOCK_DATES)
    mean_spread = _mean(list(spread.values()))
    long_excess = _mean([row["long_excess"] for row in rows])
    short_excess = _mean([row["short_excess"] for row in rows])
    return {
        "G1_net_magnitude": {"pass": mean_spread >= MAGNITUDE_THRESHOLD, "mean_spread_net": mean_spread,
                             "threshold": MAGNITUDE_THRESHOLD},
        "G2_leg_contribution": {"pass": long_excess > 0.0 and short_excess > 0.0,
                                "mean_long_excess": long_excess, "mean_short_excess": short_excess,
                                "threshold": 0.0},
        "G3_quarter_stability": {
            "pass": (stability["positive_quarters"] * QUARTER_STABILITY_DENOMINATOR
                     >= stability["included_quarters"] * QUARTER_STABILITY_NUMERATOR),
            **stability, "threshold": QUARTER_STABILITY_NUMERATOR / QUARTER_STABILITY_DENOMINATOR},
        "G4_dependence_robustness": {"pass": p_value < LOOP_ALPHA, "one_sided_block_bootstrap_p_spread_net": p_value,
                                     "alpha": LOOP_ALPHA, "replicates": BOOTSTRAP_REPLICATES,
                                     "seed": BOOTSTRAP_SEED, "block_dates": BOOTSTRAP_BLOCK_DATES},
    }


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
    counters = {key: 0 for key in ("kline-history-incomplete", "below-liquidity-floor",
                                   "missing-open-interest", "archive-without-end-of-day-row")}
    rejected = {key: 0 for key in REJECTION_KEYS}
    records = parent.load_day_records(kline_inventory, start_ms - ELIGIBILITY_DAYS * DAY_MS,
                                      end_ms + 2 * HOUR_MS)
    open_interest = load_end_of_day_open_interest(inventory, universe, start_ms - 2 * DAY_MS,
                                                  end_ms, counters)
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
            {row["date"]: row["spread_gross"] for row in rows}, replicates=BOOTSTRAP_REPLICATES,
            seed=BOOTSTRAP_SEED, block_dates=BOOTSTRAP_BLOCK_DATES) if bootstrappable(rows) else None,
        "one_hour_lagged_net_spread": {
            "dates": len(lag_values), "mean": _mean(lag_values) if lag_values else None,
            "one_sided_block_bootstrap_p": block_bootstrap_p(
                {row["date"]: row["lag_spread_net"] for row in rows if row["lag_spread_net"] is not None},
                replicates=BOOTSTRAP_REPLICATES, seed=BOOTSTRAP_SEED,
                block_dates=BOOTSTRAP_BLOCK_DATES) if bootstrappable(lag_values) else None},
        "eligibility_exclusions": counters, "rejected": rejected,
        "funding_timestamp_crossings": {"per_position_upper_bound": FUNDING_CROSSINGS_PER_POSITION,
                                        "total_upper_bound": FUNDING_CROSSINGS_PER_POSITION * positions,
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
