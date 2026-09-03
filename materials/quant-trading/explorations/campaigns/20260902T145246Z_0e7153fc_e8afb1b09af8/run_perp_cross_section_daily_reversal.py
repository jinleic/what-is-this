#!/usr/bin/env python3
"""Frozen two-window test of one-day cross-sectional reversal on Binance USDT perpetuals.

All frozen hashes and upstream archive checks are verified before candidate
parsing. Results use exclusive-create semantics. A pass is replicated
net-of-taker-fee expected-return evidence on frozen historical calendars at a
bar-open execution proxy, not deployable alpha.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import random
import zipfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
EXPLORATIONS = WORKSPACE_ROOT / "quant-trading/explorations"
CONTRACT_PATH = EXPLORATIONS / "pre_statement_perp_cross_section_reversal.json"
INVENTORY_PATH = EXPLORATIONS / "perp-cross-section-inputs.json"
AUDIT_PATH = EXPLORATIONS / "perp-cross-section-prior-window-audit.json"
POLICY_PATH = EXPLORATIONS / "frontier-loop-policy.json"

HOUR_MS = 3_600_000
DAY_MS = 86_400_000
HOURS_PER_DAY = 24
FULL_DAY_MASK = (1 << HOURS_PER_DAY) - 1
TAKER_FEE = 0.0005
ELIGIBILITY_DAYS = 20
LIQUIDITY_FLOOR_USD = 5_000_000.0
MINIMUM_CROSS_SECTION = 50
DECILE_DENOMINATOR = 10
MAGNITUDE_THRESHOLD = 0.0005
QUARTER_STABILITY_NUMERATOR = 3
QUARTER_STABILITY_DENOMINATOR = 5
FUNDING_CROSSINGS_PER_POSITION = 3
MICROSECOND_THRESHOLD = 100_000_000_000_000
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_SEED = 20_260_902
BOOTSTRAP_BLOCK_DATES = 7
EXPECTED_HEADER = (
    "open_time", "open", "high", "low", "close", "volume", "close_time",
    "quote_volume", "count", "taker_buy_volume", "taker_buy_quote_volume",
    "ignore",
)
FROZEN_PARAMETERS = {
    "taker_fee_per_side": TAKER_FEE,
    "eligibility_complete_days": ELIGIBILITY_DAYS,
    "liquidity_floor_mean_quote_volume_usd": LIQUIDITY_FLOOR_USD,
    "minimum_cross_section": MINIMUM_CROSS_SECTION,
    "decile_denominator": DECILE_DENOMINATOR,
    "magnitude_threshold_per_date": MAGNITUDE_THRESHOLD,
    "quarter_stability_fraction": [QUARTER_STABILITY_NUMERATOR,
                                   QUARTER_STABILITY_DENOMINATOR],
    "funding_timestamp_crossings_per_position": FUNDING_CROSSINGS_PER_POSITION,
    "bootstrap_replicates": BOOTSTRAP_REPLICATES,
    "bootstrap_seed": BOOTSTRAP_SEED,
    "bootstrap_block_dates": BOOTSTRAP_BLOCK_DATES,
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_iso_ms(value: str) -> int:
    return int(datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()
               * 1000)


def iso_ms(ms: int) -> str:
    return datetime.fromtimestamp(ms / 1000, timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S") + "Z"


def date_label(day_ms: int) -> str:
    return datetime.fromtimestamp(day_ms / 1000, timezone.utc).date().isoformat()


def normalize_epoch(raw: str) -> int:
    value = int(raw)
    return value // 1000 if value > MICROSECOND_THRESHOLD else value


def write_json_new(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _inside_root(path: Path, root: Path) -> Path:
    resolved_root = root.resolve()
    resolved = path.resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"input escapes workspace root: {path}") from exc
    return resolved


def _first_row_form(row: list[str]) -> bool:
    """True for the frozen header row, False for a valid headerless data row."""
    if tuple(row) == EXPECTED_HEADER:
        return True
    if len(row) == len(EXPECTED_HEADER) and row[0].isdigit():
        return False
    raise ValueError("first-row form is neither the frozen header nor a data row")


def verify_archive(entry: dict, root: Path = WORKSPACE_ROOT) -> str:
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
    try:
        observed_form = _first_row_form(first)
    except ValueError as exc:
        raise ValueError(f"{exc}: {entry['path']}") from exc
    if observed_form is not bool(entry["header_present"]):
        raise ValueError(
            f"first-row form mismatch against inventory header_present: "
            f"{entry['path']}")
    if entry["header"] != ",".join(EXPECTED_HEADER):
        raise ValueError(f"inventory header mismatch: {entry['path']}")
    return members[0]


def check_contract_parameters(contract: dict) -> None:
    frozen = contract["frozen_rule"]["parameters"]
    if frozen != FROZEN_PARAMETERS:
        raise ValueError(
            f"contract parameters differ from runner constants: "
            f"contract {frozen} runner {FROZEN_PARAMETERS}")


def verify_static_inputs(contract: dict) -> tuple[dict, dict]:
    check_contract_parameters(contract)
    expected = {
        "input_inventory_sha256": contract["assets_and_data"]
        ["input_inventory_sha256"],
        "prior_window_audit_sha256": contract["novelty_and_consumption"]
        ["prior_window_audit_sha256"],
        "loop_policy_sha256": contract["loop_multiplicity"]["policy_sha256"],
    }
    observed = {
        "input_inventory_sha256": sha256_file(INVENTORY_PATH),
        "prior_window_audit_sha256": sha256_file(AUDIT_PATH),
        "loop_policy_sha256": sha256_file(POLICY_PATH),
    }
    for field, value in expected.items():
        if observed[field] != value:
            raise ValueError(
                f"{field} mismatch: expected {value}, observed {observed[field]}")

    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
    if inventory.get("candidate_outcomes_accessed") is not False:
        raise ValueError("input inventory does not certify pre-outcome construction")
    symbols = inventory.get("symbols", [])
    if not symbols or symbols != sorted(set(symbols)):
        raise ValueError("input inventory symbols must be a sorted unique list")
    entries = inventory.get("files", [])
    if inventory.get("archive_count") != len(entries) or not entries:
        raise ValueError("input inventory archive count mismatch")
    paths = [entry["path"] for entry in entries]
    if len(paths) != len(set(paths)):
        raise ValueError("input inventory contains duplicate paths")
    symbol_set = set(symbols)
    for entry in entries:
        symbol, _ = _entry_symbol_month(entry)
        if symbol not in symbol_set:
            raise ValueError(f"inventory archive outside universe: {entry['path']}")
        verify_archive(entry)
    return inventory, {**observed, "archives_verified": len(entries)}


def expected_dates(start_ms: int, end_ms: int) -> list[int]:
    if start_ms % DAY_MS or end_ms % DAY_MS:
        raise ValueError("window boundaries must be UTC midnights")
    if end_ms <= start_ms:
        raise ValueError("window end must be after start")
    return list(range(start_ms, end_ms, DAY_MS))


def quarter_label(date: str) -> str:
    year, month = date[:4], int(date[5:7])
    return f"{year}Q{(month - 1) // 3 + 1}"


def window_months(start_ms: int, end_ms: int) -> list[str]:
    cursor = datetime.fromtimestamp(start_ms / 1000, timezone.utc).replace(
        day=1, hour=0, minute=0, second=0, microsecond=0)
    out: list[str] = []
    while int(cursor.timestamp() * 1000) < end_ms:
        out.append(cursor.strftime("%Y-%m"))
        cursor = (cursor.replace(day=28) + timedelta(days=4)).replace(day=1)
    return out


def _entry_symbol_month(entry: dict) -> tuple[str, str]:
    name = Path(entry["path"]).name
    try:
        symbol, suffix = name.split("-1h-", 1)
    except ValueError as exc:
        raise ValueError(f"unexpected inventory filename: {name}") from exc
    if not suffix.endswith(".zip"):
        raise ValueError(f"unexpected inventory filename: {name}")
    return symbol, suffix[:-4]


def _new_day_record() -> dict:
    return {"hours": 0, "duplicate": False, "bad_price": False,
            "open0": None, "open1": None, "close23": None,
            "last_hour": -1, "last_close": None, "quote_volume": 0.0,
            "trades": 0}


def _parse_archive(entry: dict, lo_ms: int, hi_ms: int,
                   records: dict[tuple[str, int], dict],
                   root: Path = WORKSPACE_ROOT) -> None:
    """Fold one archive's hourly bars into per-(symbol, day) records.

    Only bars with lo_ms <= open_time < hi_ms are read. A duplicated hour keeps
    its first occurrence and flags the day; a nonpositive or non-finite open or
    close flags the day. Neither aborts the invocation.
    """
    symbol, _ = _entry_symbol_month(entry)
    path = root / entry["path"]
    with zipfile.ZipFile(path) as archive:
        with archive.open(entry["zip_member"]) as raw:
            reader = csv.reader(io.TextIOWrapper(raw, encoding="utf-8"))
            first = next(reader)
            if _first_row_form(first) is not bool(entry["header_present"]):
                raise ValueError(f"schema first-row drift: {entry['path']}")
            rows = reader if entry["header_present"] else _chain_first(first, reader)
            for line_number, row in enumerate(rows, 2):
                if len(row) != len(EXPECTED_HEADER):
                    raise ValueError(
                        f"schema column count at {entry['path']}:{line_number}")
                try:
                    open_ms = normalize_epoch(row[0])
                    open_price = float(row[1])
                    close_price = float(row[4])
                    quote_volume = float(row[7])
                    trades = int(float(row[8]))
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"schema value at {entry['path']}:{line_number}") from exc
                if not lo_ms <= open_ms < hi_ms:
                    continue
                if open_ms % HOUR_MS != 0:
                    raise ValueError(
                        f"non-hour timestamp at {entry['path']}:{line_number}")
                day_ms = open_ms - open_ms % DAY_MS
                hour = (open_ms - day_ms) // HOUR_MS
                record = records.setdefault((symbol, day_ms), _new_day_record())
                bit = 1 << hour
                if record["hours"] & bit:
                    record["duplicate"] = True
                    continue
                record["hours"] |= bit
                if not (math.isfinite(open_price) and open_price > 0.0
                        and math.isfinite(close_price) and close_price > 0.0):
                    record["bad_price"] = True
                if not math.isfinite(quote_volume) or quote_volume < 0.0:
                    record["bad_price"] = True
                if hour == 0:
                    record["open0"] = open_price
                elif hour == 1:
                    record["open1"] = open_price
                if hour == HOURS_PER_DAY - 1:
                    record["close23"] = close_price
                if hour > record["last_hour"]:
                    record["last_hour"] = hour
                    record["last_close"] = close_price
                record["quote_volume"] += quote_volume
                record["trades"] += trades


def _chain_first(first: list[str], reader):
    yield first
    yield from reader


def load_day_records(inventory: dict, lo_ms: int,
                     hi_ms: int) -> dict[tuple[str, int], dict]:
    """Fold every inventory archive whose month touches [lo_ms, hi_ms)."""
    wanted_months = set(window_months(lo_ms, hi_ms))
    seen: set[tuple[str, str]] = set()
    records: dict[tuple[str, int], dict] = {}
    for entry in inventory["files"]:
        symbol, month = _entry_symbol_month(entry)
        if month not in wanted_months:
            continue
        if (symbol, month) in seen:
            raise ValueError(f"duplicate inventory month for {symbol} {month}")
        seen.add((symbol, month))
        _parse_archive(entry, lo_ms, hi_ms, records)
    return records


def is_complete(record: dict | None) -> bool:
    return (record is not None and record["hours"] == FULL_DAY_MASK
            and not record["duplicate"] and not record["bad_price"])


def _positive(value: float | None) -> float | None:
    if value is None or not math.isfinite(value) or value <= 0.0:
        return None
    return value


def exit_price(records: dict, symbol: str, day_ms: int,
               lagged: bool = False) -> tuple[float | None, str | None]:
    """Frozen exit chain; returns (price, fallback-key or None)."""
    prefix = "lag-" if lagged else ""
    next_day = records.get((symbol, day_ms + DAY_MS))
    if lagged and next_day is not None:
        lagged_open = _positive(next_day["open1"])
        if lagged_open is not None:
            return lagged_open, None
    if next_day is not None:
        next_open = _positive(next_day["open0"])
        if next_open is not None:
            return next_open, "lag-exit-at-next-midnight" if lagged else None
    today = records.get((symbol, day_ms))
    if today is not None:
        last_close = _positive(today["last_close"])
        if last_close is not None:
            return last_close, prefix + "early-exit-last-bar"
    return None, None


def net_return(gross: float, side: str) -> float:
    if side == "long":
        return gross * (1.0 - TAKER_FEE) - 2.0 * TAKER_FEE
    if side == "short":
        return -gross * (1.0 + TAKER_FEE) - 2.0 * TAKER_FEE
    raise ValueError(f"unknown side {side!r}")


def eligible_symbols(records: dict, symbols: list[str], day_ms: int) -> dict[str, float]:
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
            continue
        if total_quote / ELIGIBILITY_DAYS < LIQUIDITY_FLOOR_USD:
            continue
        previous = records[(symbol, day_ms - DAY_MS)]
        eligible[symbol] = previous["close23"] / previous["open0"] - 1.0
    return eligible


def build_date(records: dict, symbols: list[str], day_ms: int,
               counters: dict[str, int]) -> dict | None:
    signals = eligible_symbols(records, symbols, day_ms)
    n = len(signals)
    if n < MINIMUM_CROSS_SECTION:
        counters["insufficient-cross-section"] += 1
        return None
    k = n // DECILE_DENOMINATOR
    ranked = sorted(signals, key=lambda symbol: (signals[symbol], symbol))
    legs = {"long": ranked[:k], "short": ranked[-k:]}

    universe_gross: list[float] = []
    entries: dict[str, tuple[float, float, str | None]] = {}
    for symbol in ranked:
        today = records.get((symbol, day_ms))
        entry = _positive(today["open0"]) if today is not None else None
        if entry is None:
            continue
        price, fallback = exit_price(records, symbol, day_ms)
        if price is None:
            continue
        gross = price / entry - 1.0
        universe_gross.append(gross)
        entries[symbol] = (entry, gross, fallback)

    positions: list[dict] = []
    for side, members in legs.items():
        for symbol in members:
            today = records.get((symbol, day_ms))
            if today is None or _positive(today["open0"]) is None:
                counters["missing-entry-bar"] += 1
                continue
            if symbol not in entries:
                counters["unvaluable-exit"] += 1
                continue
            entry, gross, fallback = entries[symbol]
            if fallback == "early-exit-last-bar":
                counters["early-exit-last-bar"] += 1
            lag_net = None
            lag_entry = _positive(today["open1"])
            if lag_entry is None:
                counters["lag-missing-entry-bar"] += 1
            else:
                lag_price, lag_fallback = exit_price(records, symbol, day_ms, lagged=True)
                if lag_price is None:
                    counters["lag-unvaluable-exit"] += 1
                else:
                    if lag_fallback is not None:
                        counters[lag_fallback] += 1
                    lag_net = net_return(lag_price / lag_entry - 1.0, side)
            positions.append({
                "date": date_label(day_ms), "symbol": symbol, "side": side,
                "signal": signals[symbol], "gross": gross,
                "net": net_return(gross, side), "lag_net": lag_net,
            })
    by_side = {side: [p for p in positions if p["side"] == side]
               for side in legs}
    if any(not rows for rows in by_side.values()):
        counters["empty-leg"] += 1
        return None
    long_net = _mean([p["net"] for p in by_side["long"]])
    short_net = _mean([p["net"] for p in by_side["short"]])
    long_gross = _mean([p["gross"] for p in by_side["long"]])
    short_gross = _mean([-p["gross"] for p in by_side["short"]])
    market = _mean(universe_gross)
    lag_sides = {side: [p["lag_net"] for p in rows if p["lag_net"] is not None]
                 for side, rows in by_side.items()}
    lag_spread = None
    if all(lag_sides.values()):
        lag_spread = (_mean(lag_sides["long"]) + _mean(lag_sides["short"])) / 2.0
    date = date_label(day_ms)
    return {
        "date": date,
        "quarter": quarter_label(date),
        "eligible": n,
        "leg_size": k,
        "long_positions": len(by_side["long"]),
        "short_positions": len(by_side["short"]),
        "market_gross": market,
        "long_net": long_net,
        "short_net": short_net,
        "spread_net": (long_net + short_net) / 2.0,
        "spread_gross": (long_gross + short_gross) / 2.0,
        "long_gross": long_gross,
        "short_gross": short_gross,
        "long_excess": long_net - market,
        "short_excess": short_net + market,
        "lag_spread_net": lag_spread,
        "positions": positions,
    }


def _draw_with_rng(dates: list[str], rng: random.Random,
                   block_dates: int) -> list[str]:
    picked: list[str] = []
    while len(picked) < len(dates):
        start = rng.randrange(len(dates))
        for offset in range(block_dates):
            picked.append(dates[(start + offset) % len(dates)])
    return picked[:len(dates)]


def block_bootstrap_p(values_by_date: dict[str, float],
                      replicates: int = BOOTSTRAP_REPLICATES,
                      seed: int = BOOTSTRAP_SEED,
                      block_dates: int = BOOTSTRAP_BLOCK_DATES) -> float:
    """One-sided circular block bootstrap p for mean <= 0 over dates."""
    dates = sorted(values_by_date)
    if len(dates) < block_dates * 3:
        raise ValueError("too few dates for frozen block bootstrap")
    rng = random.Random(seed)
    at_or_below_zero = 0
    count = len(dates)
    for _ in range(replicates):
        picked = _draw_with_rng(dates, rng, block_dates)
        if sum(values_by_date[date] for date in picked) / count <= 0.0:
            at_or_below_zero += 1
    return (1 + at_or_below_zero) / (replicates + 1)


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _quarter_stability(rows: list[dict], field: str) -> dict:
    by_quarter: dict[str, list[float]] = {}
    for row in rows:
        by_quarter.setdefault(row["quarter"], []).append(row[field])
    means = {quarter: _mean(values)
             for quarter, values in sorted(by_quarter.items())}
    positive = sum(value > 0.0 for value in means.values())
    return {
        "quarterly_means": means,
        "included_quarters": len(means),
        "positive_quarters": positive,
        "positive_fraction": positive / len(means) if means else None,
    }


def run_window(contract: dict, label: str, inventory: dict) -> dict:
    window = contract["windows"][label]
    minimum = contract["minimum_sample"][label]
    start_ms = parse_iso_ms(window["start_utc_inclusive"])
    end_ms = parse_iso_ms(window["end_utc_exclusive"])
    dates = expected_dates(start_ms, end_ms)
    symbols = list(inventory["symbols"])
    lo_ms = start_ms - ELIGIBILITY_DAYS * DAY_MS
    hi_ms = end_ms + 2 * HOUR_MS
    records = load_day_records(inventory, lo_ms, hi_ms)

    counters = {key: 0 for key in (
        "insufficient-cross-section", "empty-leg", "missing-entry-bar",
        "unvaluable-exit", "early-exit-last-bar", "lag-missing-entry-bar",
        "lag-unvaluable-exit", "lag-exit-at-next-midnight",
        "lag-early-exit-last-bar")}
    day_rows: list[dict] = []
    for day_ms in dates:
        row = build_date(records, symbols, day_ms, counters)
        if row is not None:
            day_rows.append(row)

    flagged_symbol_days = {
        "duplicate-timestamp-bar": sum(
            1 for (symbol, day_ms), record in records.items()
            if record["duplicate"] and start_ms <= day_ms < end_ms),
        "nonpositive-or-nonfinite-bar": sum(
            1 for (symbol, day_ms), record in records.items()
            if record["bad_price"] and start_ms <= day_ms < end_ms),
    }
    positions = [p for row in day_rows for p in row["positions"]]
    included_quarters = len({row["quarter"] for row in day_rows})
    sample = {
        "expected_dates": len(dates),
        "included_dates": len(day_rows),
        "positions": len(positions),
        "positions_by_side": {
            side: sum(p["side"] == side for p in positions)
            for side in ("long", "short")},
        "included_quarters": included_quarters,
        "mean_eligible": _mean([row["eligible"] for row in day_rows]) if day_rows else None,
        "mean_leg_size": _mean([row["leg_size"] for row in day_rows]) if day_rows else None,
        "distinct_symbols_held": len({p["symbol"] for p in positions}),
        "rejected": counters,
        "flagged_symbol_days_in_window": flagged_symbol_days,
        "minimum": minimum,
    }
    sample_ok = (
        len(day_rows) >= minimum["included_dates"]
        and len(positions) >= minimum["positions"]
        and included_quarters >= minimum["included_quarters"]
    )
    if not sample_ok:
        return {
            "window": window,
            "status": "inconclusive",
            "reason": "frozen minimum sample not met",
            "sample": sample,
            "gate_results": None,
        }

    spread_by_date = {row["date"]: row["spread_net"] for row in day_rows}
    gross_by_date = {row["date"]: row["spread_gross"] for row in day_rows}
    lag_by_date = {row["date"]: row["lag_spread_net"] for row in day_rows
                   if row["lag_spread_net"] is not None}
    pooled = {
        "mean_spread_net": _mean(list(spread_by_date.values())),
        "mean_spread_gross": _mean(list(gross_by_date.values())),
        "mean_long_net": _mean([row["long_net"] for row in day_rows]),
        "mean_short_net": _mean([row["short_net"] for row in day_rows]),
        "mean_long_excess": _mean([row["long_excess"] for row in day_rows]),
        "mean_short_excess": _mean([row["short_excess"] for row in day_rows]),
        "mean_market_gross": _mean([row["market_gross"] for row in day_rows]),
        "spread_net_sd": _sd(list(spread_by_date.values())),
    }
    stability = _quarter_stability(day_rows, "spread_net")
    p_spread = block_bootstrap_p(spread_by_date)
    alpha = contract["loop_multiplicity"]["alpha_k"]

    gates = {
        "G1_net_magnitude": {
            "mean_spread_net": pooled["mean_spread_net"],
            "threshold": MAGNITUDE_THRESHOLD,
            "pass": pooled["mean_spread_net"] >= MAGNITUDE_THRESHOLD,
        },
        "G2_leg_contribution": {
            "mean_long_excess": pooled["mean_long_excess"],
            "mean_short_excess": pooled["mean_short_excess"],
            "threshold": 0.0,
            "pass": (pooled["mean_long_excess"] > 0.0
                     and pooled["mean_short_excess"] > 0.0),
        },
        "G3_quarter_stability": {
            "spread_net": stability,
            "threshold": QUARTER_STABILITY_NUMERATOR / QUARTER_STABILITY_DENOMINATOR,
            # Exact rational comparison: positive/included >= 3/5.
            "pass": (stability["positive_quarters"] * QUARTER_STABILITY_DENOMINATOR
                     >= stability["included_quarters"] * QUARTER_STABILITY_NUMERATOR),
        },
        "G4_dependence_robustness": {
            "one_sided_block_bootstrap_p_spread_net": p_spread,
            "alpha": alpha,
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": BOOTSTRAP_SEED,
            "block_dates": BOOTSTRAP_BLOCK_DATES,
            "pass": p_spread < alpha,
        },
    }
    all_pass = all(gate["pass"] for gate in gates.values())
    diagnostics = {
        "non_claim": contract["diagnostics"]["non_claim"],
        "gross_spread": {
            "mean": pooled["mean_spread_gross"],
            "one_sided_block_bootstrap_p": block_bootstrap_p(gross_by_date),
            "quarterly_means": _quarter_stability(day_rows, "spread_gross")["quarterly_means"],
        },
        "one_hour_lagged_net_spread": {
            "dates": len(lag_by_date),
            "mean": _mean(list(lag_by_date.values())) if lag_by_date else None,
            "one_sided_block_bootstrap_p": (
                block_bootstrap_p(lag_by_date)
                if len(lag_by_date) >= BOOTSTRAP_BLOCK_DATES * 3 else None),
        },
        "funding_timestamp_crossings": {
            "per_position_upper_bound": FUNDING_CROSSINGS_PER_POSITION,
            "total_upper_bound": FUNDING_CROSSINGS_PER_POSITION * len(positions),
            "subtracted": False,
        },
        "yearly_mean_spread_net": _yearly_means(day_rows, "spread_net"),
        "mean_leg_gross": {
            "long": _mean([row["long_gross"] for row in day_rows]),
            "short": _mean([row["short_gross"] for row in day_rows]),
        },
    }
    return {
        "window": window,
        "status": "pass-on-this-sample" if all_pass else "falsified",
        "reason": None if all_pass else "at least one frozen gate is false",
        "sample": sample,
        "pooled": pooled,
        "cost_model": {"taker_fee_per_side": TAKER_FEE,
                       "round_trip_fee_drag": 2 * TAKER_FEE,
                       "funding_subtracted": False,
                       "slippage_or_impact": None},
        "gate_results": gates,
        "gate_descriptions": contract["gates"],
        "diagnostics": diagnostics,
        "daily": [{key: value for key, value in row.items() if key != "positions"}
                  for row in day_rows],
    }


def _sd(values: list[float]) -> float | None:
    if len(values) < 2:
        return None
    mean = _mean(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))


def _yearly_means(rows: list[dict], field: str) -> dict[str, float]:
    by_year: dict[str, list[float]] = {}
    for row in rows:
        by_year.setdefault(row["date"][:4], []).append(row[field])
    return {year: _mean(values) for year, values in sorted(by_year.items())}


def validate_primary_artifact(primary: dict, expected: dict) -> None:
    if primary.get("role") != "primary":
        raise ValueError("replication artifact is not role primary")
    for field, value in expected.items():
        if primary.get(field) != value:
            raise ValueError(
                f"replication primary {field} mismatch: expected {value}, "
                f"observed {primary.get(field)}")
    status = primary.get("result", {}).get("status")
    if status != "pass-on-this-sample":
        raise ValueError(
            "replication requires a committed primary pass; "
            f"observed {status!r}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT_PATH)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--window", choices=("primary", "replication"),
                        required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--primary-artifact", type=Path)
    parser.add_argument("--primary-artifact-sha256")
    parser.add_argument("--verify-inputs-only", action="store_true")
    args = parser.parse_args(argv)

    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    contract_sha = sha256_file(args.contract.resolve())
    if contract_sha != args.contract_sha256:
        raise SystemExit(
            f"contract hash mismatch: expected {args.contract_sha256}, "
            f"observed {contract_sha}")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    try:
        inventory, integrity = verify_static_inputs(contract)
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise SystemExit(f"static integrity failure: {exc}") from exc

    runner_sha = sha256_file(Path(__file__).resolve())
    expected_primary = {
        "contract_sha256": contract_sha,
        "input_inventory_sha256": integrity["input_inventory_sha256"],
        "prior_window_audit_sha256": integrity["prior_window_audit_sha256"],
        "loop_policy_sha256": integrity["loop_policy_sha256"],
        "runner_sha256": runner_sha,
    }
    if args.verify_inputs_only:
        print(json.dumps({**expected_primary,
                          "archives_verified": integrity["archives_verified"]},
                         indent=2, sort_keys=True))
        return 0

    primary_sha = None
    if args.window == "replication":
        if args.primary_artifact is None or args.primary_artifact_sha256 is None:
            raise SystemExit(
                "replication requires --primary-artifact and "
                "--primary-artifact-sha256")
        if not args.primary_artifact.is_file():
            raise SystemExit(
                f"replication primary artifact absent: {args.primary_artifact}")
        primary_sha = sha256_file(args.primary_artifact)
        if primary_sha != args.primary_artifact_sha256:
            raise SystemExit(
                "replication primary artifact hash does not match committed hash")
        primary = json.loads(args.primary_artifact.read_text(encoding="utf-8"))
        try:
            validate_primary_artifact(primary, expected_primary)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc

    try:
        result = run_window(contract, args.window, inventory)
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise SystemExit(f"outcome input/schema failure: {exc}") from exc
    now_ms = int(datetime.now(timezone.utc).timestamp() * 1000)
    payload = {
        "schema_version": 1,
        "id": contract["id"],
        "created_at_utc": iso_ms(now_ms),
        **expected_primary,
        "archives_verified": integrity["archives_verified"],
        "role": args.window,
        "result": result,
        "claim_boundary": contract["claim_boundary"],
    }
    if args.window == "replication":
        payload["primary_artifact"] = str(
            args.primary_artifact.resolve().relative_to(WORKSPACE_ROOT.resolve()))
        payload["primary_artifact_sha256"] = primary_sha
    write_json_new(args.output, payload)
    print(json.dumps({
        f"{args.window}_status": result["status"],
        "included_dates": result["sample"]["included_dates"],
        "positions": result["sample"]["positions"],
        "output": str(args.output),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
