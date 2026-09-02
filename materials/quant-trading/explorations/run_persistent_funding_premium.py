#!/usr/bin/env python3
"""Contract-compliant preregistered test: persistent funding premium.

Frozen contract:
    quant-trading/explorations/persistent-funding-premium-contract.json
    sha256 04319eaef33f933c767b20dd896c711b065b7d02861f6861cb0c910227d3ad5e

Rule summary (frozen before outcome access):
  * Signal: 18 consecutive strictly positive funding events (events i-17..i,
    including i). Entry at the first 1h spot bar open strictly after the
    event and within 2h of it; skip and count otherwise. Long spot 0.5x NAV
    + short perpetual 0.5x NAV; gross 1.0x; both legs pay 7 bps per side at
    entry and exit.
  * Exit at the 1h bar open 240 hours after the entry bar open (not the
    close). Funding accrues on the short perp leg at each exchange boundary
    strictly after entry and at-or-before exit; signed exactly.
  * Hard liquidation gate (frozen): if within the holding period the perp
    HIGH exceeds 2x the perp entry price (short-leg cash floor breach on
    0.5x collateral), unwind early at that hour's perp HIGH and the same
    hour's spot LOW, with the same entry/exit fee schedule.
  * Gates: minima met (>=8 pooled, >=4 per asset) else inconclusive;
    pooled mean and pooled median > 0 else falsified; worst entry <= -0.10
    rejects as unstable regardless of mean/median.
  * A pass triggers exactly one replication run on 2025-07-01..2026-01-01
    with identical gates. No third window, no parameter search, no pooling.
  * Sample disjointness (fresh data): 2024-07..2024-12 dev, 2025-01..2025-06
    holdout, 2025-07..2025-12 replication. No prior mechanism consumed
    these windows: funding-crowding used 2020-08..2024-12, carry used the
    same plus an invalid retrospective 2025 diagnostic, and all
    focused-campaign screens used live 2026-08 captures only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import zipfile
from bisect import bisect_left, bisect_right
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = WORKSPACE_ROOT / (
    "quant-trading/explorations/persistent-funding-premium-contract.json"
)
FROZEN_CONTRACT_SHA256 = (
    "04319eaef33f933c767b20dd896c711b065b7d02861f6861cb0c910227d3ad5e"
)
ASSETS = ["BTCUSDT", "ETHUSDT"]
HOLDING_HOURS = 240
ENTRY_LAG_LIMIT_MS = 2 * 60 * 60 * 1000
EXIT_LAG_LIMIT_MS = 2 * 60 * 60 * 1000
TAIL_FLOOR = -0.10
SPOT_NAV_FRACTION = 0.5
PERP_NAV_FRACTION = 0.5


# --------------------------------------------------------------------------
# utilities
# --------------------------------------------------------------------------

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iso_ms(ms: int) -> str:
    return datetime.fromtimestamp(ms / 1000, timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S.%f"
    )[:-3] + "Z"


def parse_iso_ms(value: str) -> int:
    return int(
        datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp() * 1000
    )


def normalize_epoch(raw: str) -> int:
    """Binance archives: ms before 2025, us from 2025. Return UTC ms."""
    value = int(raw)
    return value // 1000 if value > 10_000_000_000_000 else value


# --------------------------------------------------------------------------
# checksum verification (upstream .CHECKSUM binding)
# --------------------------------------------------------------------------

def verify_checksum(rel_path: str) -> str:
    path = WORKSPACE_ROOT / rel_path
    checksum_path = path.with_name(path.name + ".CHECKSUM")
    expected = checksum_path.read_text(encoding="utf-8").split()[0].lower()
    observed = sha256_file(path)
    if observed != expected:
        raise ValueError(
            f"upstream checksum mismatch for {rel_path}: {observed} != {expected}"
        )
    return observed


# --------------------------------------------------------------------------
# frozen monthly archive readers
# --------------------------------------------------------------------------

def month_numbers(window: dict) -> list[tuple[int, int]]:
    """Months intersecting [start, end_exclusive), plus one context month
    before start for funding signals (never for bars inside gates)."""
    start = parse_iso_ms(window["start"])
    end = parse_iso_ms(window["end_exclusive"])
    months: set[tuple[int, int]] = set()

    def add(ms: int) -> None:
        dt = datetime.fromtimestamp(ms / 1000, timezone.utc)
        months.add((dt.year, dt.month))

    # one context month before start (signal history only)
    month_start = datetime.fromtimestamp(start / 1000, timezone.utc).replace(
        day=1, hour=0, minute=0, second=0, microsecond=0
    )
    add(int(month_start.timestamp() * 1000) - 1)
    cursor = start
    while cursor < end:
        add(cursor)
        # advance one month
        dt = datetime.fromtimestamp(cursor / 1000, timezone.utc)
        if dt.month == 12:
            dt = dt.replace(year=dt.year + 1, month=1)
        else:
            dt = dt.replace(month=dt.month + 1)
        cursor = int(dt.replace(day=1, hour=0, minute=0, second=0, microsecond=0).timestamp() * 1000)
    return sorted(months)


def read_zip_csv(rel_path: str) -> list[list[str]]:
    path = WORKSPACE_ROOT / rel_path
    with zipfile.ZipFile(path) as archive:
        name = archive.namelist()[0]
        with archive.open(name) as handle:
            text = io.TextIOWrapper(handle, encoding="utf-8").read()
    return list(csv.reader(io.StringIO(text)))


def load_funding(asset: str, window: dict) -> tuple[list[int], list[float]]:
    """Funding events sorted by exchange calc_time. Includes one context
    month before the window for trailing-signal history only."""
    times: list[tuple[int, int]] = []
    rates: list[tuple[int, float]] = []
    for year, month in month_numbers(window):
        rel = (
            f"quant-trading/data/raw/binance/funding-rate/{asset}/"
            f"{asset}-fundingRate-{year:04d}-{month:02d}.zip"
        )
        if not (WORKSPACE_ROOT / rel).exists():
            continue
        verify_checksum(rel)
        for index, row in enumerate(read_zip_csv(rel)[1:]):
            calc_time = int(row[0])
            times.append((calc_time, index))
            rates.append((calc_time, float(row[2])))
    times.sort()
    rates.sort()
    return [t for t, _ in times], [r for _, r in rates]


def load_bars(asset: str, kind: str, window: dict) -> dict:
    """1h bars sorted by exchange open_time: times, open, high, low.
    Bars are full archive rows; gate math filters by window boundaries."""
    rows: list[tuple[int, int]] = []
    opens: list[tuple[int, float]] = []
    highs: list[tuple[int, float]] = []
    lows: list[tuple[int, float]] = []
    start_ms = parse_iso_ms(window["start"])
    end_ms = parse_iso_ms(window["end_exclusive"])
    for year, month in month_numbers(window):
        rel = (
            f"quant-trading/data/raw/binance/{kind}/{asset}/"
            f"{asset}-1h-{year:04d}-{month:02d}.zip"
        )
        if not (WORKSPACE_ROOT / rel).exists():
            continue
        verify_checksum(rel)
        for row in read_zip_csv(rel)[1:]:
            open_ms = normalize_epoch(row[0])
            # keep bars overlapping the window plus one hour of context
            if open_ms > end_ms + 3600 * 1000 or open_ms < start_ms - 49 * 3600 * 1000:
                continue
            rows.append((open_ms, 0))
            opens.append((open_ms, float(row[1])))
            highs.append((open_ms, float(row[2])))
            lows.append((open_ms, float(row[3])))
    rows.sort()
    opens.sort(key=lambda item: item[0])
    highs.sort(key=lambda item: item[0])
    lows.sort(key=lambda item: item[0])
    return {
        "times": [t for t, _ in rows],
        "open": [v for _, v in opens],
        "high": [v for _, v in highs],
        "low": [v for _, v in lows],
    }


# --------------------------------------------------------------------------
# window data (loaded once per window)
# --------------------------------------------------------------------------

class WindowData:
    def __init__(self, window: dict, label: str) -> None:
        self.window = window
        self.label = label
        self.funding: dict[str, tuple[list[int], list[float]]] = {}
        self.spot: dict[str, dict] = {}
        self.perp: dict[str, dict] = {}
        for asset in ASSETS:
            self.funding[asset] = load_funding(asset, window)
            self.spot[asset] = load_bars(asset, "spot-1h", window)
            self.perp[asset] = load_bars(asset, "um-perp-1h", window)


# --------------------------------------------------------------------------
# frozen trade engine (one persistence variant over one window)
# --------------------------------------------------------------------------

def trailing_all_positive(rates: list[float], index: int, count: int) -> bool:
    if index + 1 < count:
        return False
    return all(value > 0 for value in rates[index + 1 - count : index + 1])


def run_variant(asset: str, data: WindowData, persistence_count: int, fee_bps: float) -> dict:
    funding_times, funding_rates = data.funding[asset]
    spot = data.spot[asset]
    perp = data.perp[asset]
    if not (spot["times"] and perp["times"]):
        raise ValueError(f"unavailable bar sets for {asset} in {data.label}")

    start_ms = parse_iso_ms(data.window["start"])
    end_ms = parse_iso_ms(data.window["end_exclusive"])
    holdout_clone = None  # unused placeholder
    fee = fee_bps / 10_000.0  # per side per leg
    cost_per_cycle = 2 * fee * (SPOT_NAV_FRACTION + PERP_NAV_FRACTION) * 2 / 2
    # entry fees on both legs + exit fees on both legs:
    cost_per_cycle = fee * (SPOT_NAV_FRACTION + PERP_NAV_FRACTION) * 2

    entries = 0
    skipped_no_bar = 0
    skipped_outside_window = 0
    forced_unwinds = 0
    detail: list[dict] = []
    position_open_until = -1

    for index, event_time in enumerate(funding_times):
        if event_time < position_open_until:
            continue
        if not trailing_all_positive(funding_rates, index, persistence_count):
            continue
        # entry bar: first spot open strictly after the event
        entry_idx = bisect_left(spot["times"], event_time + 1)
        if entry_idx >= len(spot["times"]):
            skipped_no_bar += 1
            continue
        entry_spot_time = spot["times"][entry_idx]
        lag = entry_spot_time - event_time
        if lag > ENTRY_LAG_LIMIT_MS:
            skipped_no_bar += 1
            continue
        if entry_spot_time < start_ms or entry_spot_time >= end_ms:
            skipped_outside_window += 1
            continue
        pidx = bisect_left(perp["times"], entry_spot_time)
        if pidx >= len(perp["times"]) or perp["times"][pidx] != entry_spot_time:
            skipped_no_bar += 1
            continue
        spot_entry = spot["open"][entry_idx]
        perp_entry = perp["open"][pidx]

        exit_target = entry_spot_time + HOLDING_HOURS * 3600 * 1000
        exit_idx = bisect_left(spot["times"], exit_target)
        if exit_idx >= len(spot["times"]):
            skipped_outside_window += 1  # outcome unavailable inside frozen sample
            continue
        exit_spot_time = spot["times"][exit_idx]
        if exit_spot_time >= end_ms:
            skipped_outside_window += 1
            continue
        if abs(exit_spot_time - exit_target) > EXIT_LAG_LIMIT_MS:
            skipped_no_bar += 1
            continue
        pexit_idx = bisect_left(perp["times"], exit_spot_time)
        if pexit_idx >= len(perp["times"]) or perp["times"][pexit_idx] != exit_spot_time:
            skipped_no_bar += 1
            continue

        # hour-by-hour hard liquidation audit until exit (frozen rule)
        forced = False
        forced_at = None
        audit_idx = bisect_right(perp["times"], entry_spot_time)
        while audit_idx < len(perp["times"]):
            t = perp["times"][audit_idx]
            if t > exit_spot_time:
                break
            if perp["high"][audit_idx] > 2.0 * perp_entry:
                forced = True
                forced_at = t
                break
            audit_idx += 1

        if forced:
            # unwind at the same-hour adverse extremes: perp HIGH, spot LOW
            sidx = bisect_left(spot["times"], forced_at)
            if sidx >= len(spot["times"]) or spot["times"][sidx] != forced_at:
                skipped_no_bar += 1
                continue
            spot_exit = spot["low"][sidx]
            perp_exit = perp["high"][audit_idx]
            exit_time = forced_at
            forced_unwinds += 1
            position_open_until = exit_time + 1
        else:
            spot_exit = spot["open"][exit_idx]
            perp_exit = perp["open"][pexit_idx]
            exit_time = exit_spot_time
            position_open_until = exit_spot_time + 1

        # price legs with frozen 0.5/0.5 NAV fractions
        spot_leg = spot_exit / spot_entry - 1.0
        perp_leg = -(perp_exit / perp_entry - 1.0)  # short perp
        price_component = SPOT_NAV_FRACTION * spot_leg + PERP_NAV_FRACTION * perp_leg

        # signed funding on the short perp leg: events strictly after entry,
        # at-or-before the realized exit time
        funding_component = 0.0
        j = bisect_right(funding_times, entry_spot_time)
        while j < len(funding_times) and funding_times[j] <= exit_time:
            funding_component += PERP_NAV_FRACTION * funding_rates[j]
            j += 1

        net_fraction = price_component + funding_component - cost_per_cycle
        entries += 1
        detail.append(
            {
                "asset": asset,
                "signal_event_ms": event_time,
                "signal_event_utc": iso_ms(event_time),
                "exit_time_ms": exit_time,
                "exit_time_utc": iso_ms(exit_time),
                "forced_unwind": forced,
                "spot_entry": spot_entry,
                "spot_exit": spot_exit,
                "perp_entry": perp_entry,
                "perp_exit": perp_exit,
                "price_component": price_component,
                "funding_component": funding_component,
                "cost_per_cycle": cost_per_cycle,
                "raw_fraction": price_component + funding_component,
                "net_fraction": net_fraction,
            }
        )

    return {
        "asset": asset,
        "persistence_count": persistence_count,
        "entries": entries,
        "skipped_no_bar": skipped_no_bar,
        "skipped_outside_window": skipped_outside_window,
        "forced_unwinds": forced_unwinds,
        "detail": detail,
    }


# --------------------------------------------------------------------------
# frozen gates
# --------------------------------------------------------------------------

def evaluate_gates(detail: list[dict], minimum_pooled: int, minimum_per_asset: int) -> dict:
    pooled = [row["net_fraction"] for row in detail]
    per_asset_counts = {
        asset: sum(1 for row in detail if row["asset"] == asset) for asset in ASSETS
    }
    minima_met = (
        len(pooled) >= minimum_pooled
        and all(c >= minimum_per_asset for c in per_asset_counts.values())
    )
    pooled_mean = sum(pooled) / len(pooled) if pooled else None
    pooled_median = median(pooled) if pooled else None
    worst = min(pooled) if pooled else None
    falsified = minima_met and (
        pooled_mean <= 0
        or pooled_median <= 0
        or worst is None
        or worst <= TAIL_FLOOR
    )
    if not minima_met:
        status = "inconclusive"
    elif falsified:
        status = "falsified"
    else:
        status = "validation-pass-on-this-sample"
    return {
        "status": status,
        "minima_met": minima_met,
        "pooled_entries": len(pooled),
        "per_asset_entries": per_asset_counts,
        "pooled_mean": pooled_mean,
        "pooled_median": pooled_median,
        "worst_entry": worst,
        "mean_gate": pooled_mean is not None and pooled_mean > 0,
        "median_gate": pooled_median is not None and pooled_median > 0,
        "tail_gate": worst is not None and worst > TAIL_FLOOR,
        "tail_floor": TAIL_FLOOR,
    }


def variant_bundle(data: WindowData, persistence_count: int, fee_bps: float, minimum_pooled: int, minimum_per_asset: int) -> dict:
    per_asset = [run_variant(asset, data, persistence_count, fee_bps) for asset in ASSETS]
    detail: list[dict] = []
    for item in per_asset:
        detail.extend(item["detail"])
    return {
        "persistence_count": persistence_count,
        "per_asset": {
            item["asset"]: {
                "entries": item["entries"],
                "skipped_no_bar": item["skipped_no_bar"],
                "skipped_outside_window": item["skipped_outside_window"],
                "forced_unwinds": item["forced_unwinds"],
            }
            for item in per_asset
        },
        "gates": evaluate_gates(detail, minimum_pooled, minimum_per_asset),
        "entries_detail": detail,
    }


# --------------------------------------------------------------------------
# orchestration: one holdout resolution, conditional single replication
# --------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-checksums", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=WORKSPACE_ROOT
        / "quant-trading/explorations/persistent-funding-premium/results.json",
    )
    parser.add_argument(
        "--replication-output",
        type=Path,
        default=WORKSPACE_ROOT
        / "quant-trading/explorations/persistent-funding-premium/replication_results.json",
    )
    args = parser.parse_args(argv)

    observed_sha = sha256_file(CONTRACT_PATH)
    if observed_sha != FROZEN_CONTRACT_SHA256:
        raise SystemExit(
            "frozen contract changed: "
            f"expected {FROZEN_CONTRACT_SHA256}, observed {observed_sha}"
        )
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    fee_bps = contract["cost_model"]["fee_bps_per_side_per_leg"]["spot"]
    primary_gate = contract["primary_falsification_gates"][0]
    minimum_pooled = primary_gate["minimum_pooled_entries"]
    minimum_per_asset = primary_gate["minimum_per_asset_entries"]
    variant_primary = contract["preregistered_variants"][0]["persistence_count"]
    variant_robustness = contract["preregistered_variants"][1]["persistence_count"]

    if args.verify_checksums:
        checked = []
        for window in (
            contract["sample"]["development_window_utc"],
            contract["sample"]["holdout_window_utc"],
            contract["replication_contract"]["second_window_utc"],
        ):
            for asset in ASSETS:
                for year, month in month_numbers(window):
                    for kind, template in (
                        ("funding", "{a}-fundingRate-{ym}.zip"),
                        ("spot-1h", "{a}-1h-{ym}.zip"),
                        ("um-perp-1h", "{a}-1h-{ym}.zip"),
                    ):
                        rel = (
                            f"quant-trading/data/raw/binance/{kind}/{asset}/"
                            + template.format(a=asset, ym=f"{year:04d}-{month:02d}")
                        )
                        if (WORKSPACE_ROOT / rel).exists():
                            verify_checksum(rel)
                            checked.append(rel)
        print(json.dumps({"verify_checksums": "ok", "files": len(checked)}, indent=2))
        return 0

    # Dev window: feasibility and direction context. If the frozen 18/18
    # rule yields zero dev entries, stop without touching holdout data.
    dev_window = contract["sample"]["development_window_utc"]
    dev_data = WindowData(dev_window, "development")
    dev_primary = variant_bundle(
        dev_data, variant_primary, fee_bps, minimum_pooled, minimum_per_asset
    )
    dev_robustness = variant_bundle(
        dev_data, variant_robustness, fee_bps, 1, 0
    )
    del dev_data

    if dev_primary["gates"]["pooled_entries"] == 0:
        payload = {
            "schema_version": 1,
            "campaign": contract["campaign"],
            "created_at_utc": iso_ms(int(datetime.now(timezone.utc).timestamp() * 1000)),
            "contract_sha256": observed_sha,
            "status": "blocked-no-dev-sample",
            "development": {
                "primary": {k: v for k, v in dev_primary.items() if k != "entries_detail"},
                "robustness_diagnostic": {k: v for k, v in dev_robustness.items() if k != "entries_detail"},
            },
            "claim_boundary": contract["claim_boundary"],
        }
        if args.output.exists():
            raise FileExistsError(f"refusing to overwrite frozen artifact: {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({"status": "blocked-no-dev-sample", "output": str(args.output)}, indent=2))
        return 0

    # Holdout resolution (single). Entries strictly inside the holdout window.
    holdout_window = contract["sample"]["holdout_window_utc"]
    holdout_data = WindowData(holdout_window, "holdout")
    holdout_primary = variant_bundle(
        holdout_data, variant_primary, fee_bps, minimum_pooled, minimum_per_asset
    )
    holdout_robustness = variant_bundle(
        holdout_data, variant_robustness, fee_bps, minimum_pooled, minimum_per_asset
    )
    del holdout_data
    primary_status = holdout_primary["gates"]["status"]

    payload = {
        "schema_version": 1,
        "campaign": contract["campaign"],
        "created_at_utc": iso_ms(int(datetime.now(timezone.utc).timestamp() * 1000)),
        "contract_sha256": observed_sha,
        "classification": contract["classification"],
        "claim_boundary": contract["claim_boundary"],
        "windows": {
            "development": dev_window,
            "holdout": holdout_window,
        },
        "development": {
            "primary": {k: v for k, v in dev_primary.items() if k != "entries_detail"},
            "robustness_diagnostic": {k: v for k, v in dev_robustness.items() if k != "entries_detail"},
        },
        "holdout": {
            "primary": {k: v for k, v in holdout_primary.items() if k != "entries_detail"},
            "robustness_diagnostic": {k: v for k, v in holdout_robustness.items() if k != "entries_detail"},
        },
        "primary_status": primary_status,
    }
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite frozen artifact: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"primary_status": primary_status, "output": str(args.output)}, indent=2))

    # Single preregistered replication, only on a primary pass, only if the
    # second-window data exists; otherwise blocked-no-second-sample-data.
    if primary_status != "validation-pass-on-this-sample":
        return 0

    second_window = contract["replication_contract"]["second_window_utc"]
    required_months = month_numbers(second_window)
    missing: list[str] = []
    for asset in ASSETS:
        for year, month in required_months:
            for kind, template in (
                ("funding", "{a}-fundingRate-{ym}.zip"),
                ("spot-1h", "{a}-1h-{ym}.zip"),
                ("um-perp-1h", "{a}-1h-{ym}.zip"),
            ):
                rel = (
                    f"quant-trading/data/raw/binance/{kind}/{asset}/"
                    + template.format(a=asset, ym=f"{year:04d}-{month:02d}")
                )
                if not (WORKSPACE_ROOT / rel).exists():
                    missing.append(rel)
    if missing:
        replication_payload = {
            "schema_version": 1,
            "campaign": contract["campaign"],
            "created_at_utc": iso_ms(int(datetime.now(timezone.utc).timestamp() * 1000)),
            "contract_sha256": observed_sha,
            "status": "blocked-no-second-sample-data",
            "missing_files": missing,
            "claim_boundary": contract["claim_boundary"],
        }
        if args.replication_output.exists():
            raise FileExistsError(f"refusing to overwrite frozen artifact: {args.replication_output}")
        args.replication_output.write_text(
            json.dumps(replication_payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({"replication_status": "blocked-no-second-sample-data", "output": str(args.replication_output)}, indent=2))
        return 0

    replication_data = WindowData(second_window, "replication")
    replication_primary = variant_bundle(
        replication_data, variant_primary, fee_bps, minimum_pooled, minimum_per_asset
    )
    replication_robustness = variant_bundle(
        replication_data, variant_robustness, fee_bps, minimum_pooled, minimum_per_asset
    )
    del replication_data
    replication_status = replication_primary["gates"]["status"]
    replication_payload = {
        "schema_version": 1,
        "campaign": contract["campaign"],
        "created_at_utc": iso_ms(int(datetime.now(timezone.utc).timestamp() * 1000)),
        "contract_sha256": observed_sha,
        "second_window_utc": second_window,
        "replication": {
            "primary": {k: v for k, v in replication_primary.items() if k != "entries_detail"},
            "robustness_diagnostic": {k: v for k, v in replication_robustness.items() if k != "entries_detail"},
        },
        "replication_status": replication_status,
        "claim_boundary": contract["claim_boundary"],
    }
    if args.replication_output.exists():
        raise FileExistsError(f"refusing to overwrite frozen artifact: {args.replication_output}")
    args.replication_output.write_text(
        json.dumps(replication_payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"replication_status": replication_status, "output": str(args.replication_output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
