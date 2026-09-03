#!/usr/bin/env python3
"""Frozen two-window test of the 21:00-23:00 UTC net-of-taker-fee premium.

All frozen hashes and upstream archive checks are verified before candidate
parsing. Results use exclusive-create semantics. A pass is replicated
net-of-cost expected-return evidence on frozen historical calendars, not
deployable alpha.
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
CONTRACT_PATH = EXPLORATIONS / "pre_statement_utc_2100_2300.json"
INVENTORY_PATH = EXPLORATIONS / "utc-2100-2300-inputs.json"
AUDIT_PATH = EXPLORATIONS / "utc-2100-2300-prior-window-audit.json"
POLICY_PATH = EXPLORATIONS / "frontier-loop-policy.json"

ASSETS = ("BTCUSDT", "ETHUSDT")
HOUR_MS = 3_600_000
DAY_MS = 86_400_000
HOURS_PER_DAY = 24
ENTRY_HOUR = 21
EXIT_HOUR = 23
WINDOW_HOURS = EXIT_HOUR - ENTRY_HOUR
TAKER_FEE = 0.0005
MICROSECOND_THRESHOLD = 100_000_000_000_000
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_SEED = 20_260_902
BOOTSTRAP_BLOCK_DATES = 7
EXPECTED_HEADER = (
    "open_time", "open", "high", "low", "close", "volume", "close_time",
    "quote_volume", "count", "taker_buy_volume", "taker_buy_quote_volume",
    "ignore",
)


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


def verify_static_inputs(contract: dict) -> tuple[dict, dict]:
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
    if tuple(inventory.get("assets", ())) != ASSETS:
        raise ValueError("input inventory asset order mismatch")
    entries = inventory.get("files", [])
    if inventory.get("archive_count") != len(entries) or not entries:
        raise ValueError("input inventory archive count mismatch")
    paths = [entry["path"] for entry in entries]
    if len(paths) != len(set(paths)):
        raise ValueError("input inventory contains duplicate paths")
    for entry in entries:
        verify_archive(entry)
    return inventory, {**observed, "archives_verified": len(entries)}


def expected_dates(start_ms: int, end_ms: int) -> list[int]:
    if start_ms % DAY_MS or end_ms % DAY_MS:
        raise ValueError("window boundaries must be UTC midnights")
    return list(range(start_ms, end_ms, DAY_MS))


def quarter_label(date: str) -> str:
    year, month = date[:4], int(date[5:7])
    return f"{year}Q{(month - 1) // 3 + 1}"


def daily_observation(asset: str, day_ms: int, bars: dict[int, dict]) -> dict:
    if len(bars) != HOURS_PER_DAY or any(
            day_ms + hour * HOUR_MS not in bars for hour in range(HOURS_PER_DAY)):
        raise ValueError("complete hourly grid missing")
    for bar in bars.values():
        for field in ("open", "close"):
            value = bar[field]
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError("nonpositive price in required bar")
    p_in = bars[day_ms + ENTRY_HOUR * HOUR_MS]["open"]
    p_out = bars[day_ms + EXIT_HOUR * HOUR_MS]["open"]
    p_0 = bars[day_ms]["open"]
    p_24 = bars[day_ms + (HOURS_PER_DAY - 1) * HOUR_MS]["close"]
    date = datetime.fromtimestamp(day_ms / 1000, timezone.utc).date().isoformat()
    gross = p_out / p_in
    return {
        "asset": asset,
        "date": date,
        "quarter": quarter_label(date),
        "net": gross * (1.0 - TAKER_FEE) - (1.0 + TAKER_FEE),
        "x": math.log(gross) - (WINDOW_HOURS / HOURS_PER_DAY) * math.log(p_24 / p_0),
    }


def window_months(start_ms: int, end_ms: int) -> list[str]:
    cursor = datetime.fromtimestamp(start_ms / 1000, timezone.utc).replace(
        day=1, hour=0, minute=0, second=0, microsecond=0)
    out: list[str] = []
    while int(cursor.timestamp() * 1000) < end_ms:
        out.append(cursor.strftime("%Y-%m"))
        cursor = (cursor.replace(day=28) + timedelta(days=4)).replace(day=1)
    return out


def _entry_asset_month(entry: dict) -> tuple[str, str]:
    name = Path(entry["path"]).name
    try:
        asset, suffix = name.split("-1h-", 1)
    except ValueError as exc:
        raise ValueError(f"unexpected inventory filename: {name}") from exc
    if not suffix.endswith(".zip"):
        raise ValueError(f"unexpected inventory filename: {name}")
    return asset, suffix[:-4]


def _parse_archive(entry: dict, start_ms: int, end_ms: int,
                   target_days: set[int], seen: set[int],
                   days: dict[int, dict[int, dict]],
                   duplicated_days: set[int],
                   root: Path = WORKSPACE_ROOT) -> None:
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
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"schema value at {entry['path']}:{line_number}") from exc
                if not start_ms <= open_ms < end_ms:
                    continue
                if open_ms % HOUR_MS != 0:
                    raise ValueError(
                        f"non-hour timestamp at {entry['path']}:{line_number}")
                day_ms = open_ms - open_ms % DAY_MS
                if open_ms in seen:
                    # Frozen rule: a duplicate bar rejects that asset-date and
                    # is counted; it never aborts the invocation.
                    duplicated_days.add(day_ms)
                    continue
                seen.add(open_ms)
                if day_ms not in target_days:
                    continue
                days.setdefault(day_ms, {})[open_ms] = {
                    "open": open_price, "close": close_price}


def _chain_first(first: list[str], reader):
    yield first
    yield from reader


def load_asset_days(inventory: dict, asset: str, start_ms: int,
                    end_ms: int) -> tuple[dict[int, dict[int, dict]], set[int]]:
    """Return (bars by target day, days holding a duplicated timestamp)."""
    wanted_months = window_months(start_ms, end_ms)
    entries_by_month: dict[str, dict] = {}
    for entry in inventory["files"]:
        entry_asset, month = _entry_asset_month(entry)
        if entry_asset == asset and month in wanted_months:
            if month in entries_by_month:
                raise ValueError(f"duplicate inventory month for {asset} {month}")
            entries_by_month[month] = entry
    if sorted(entries_by_month) != sorted(wanted_months):
        raise ValueError(
            f"inventory months for {asset}: expected {wanted_months}, "
            f"observed {sorted(entries_by_month)}")

    target_days = set(expected_dates(start_ms, end_ms))
    days: dict[int, dict[int, dict]] = {}
    seen: set[int] = set()
    duplicated_days: set[int] = set()
    for month in wanted_months:
        _parse_archive(entries_by_month[month], start_ms, end_ms,
                       target_days, seen, days, duplicated_days)
    return days, duplicated_days


def _draw_with_rng(dates: list[str], rng: random.Random,
                   block_dates: int) -> list[str]:
    picked: list[str] = []
    while len(picked) < len(dates):
        start = rng.randrange(len(dates))
        for offset in range(block_dates):
            picked.append(dates[(start + offset) % len(dates)])
    return picked[:len(dates)]


def block_bootstrap_p(rows: list[dict], field: str,
                      replicates: int = BOOTSTRAP_REPLICATES,
                      seed: int = BOOTSTRAP_SEED,
                      block_dates: int = BOOTSTRAP_BLOCK_DATES) -> float:
    by_date: dict[str, list[float]] = {}
    for row in rows:
        by_date.setdefault(row["date"], []).append(row[field])
    dates = sorted(by_date)
    if len(dates) < block_dates * 3:
        raise ValueError("too few dates for frozen block bootstrap")
    if any(len(by_date[date]) != len(ASSETS) for date in dates):
        raise ValueError("bootstrap rows are not asset-paired by date")
    date_sums = {date: sum(by_date[date]) for date in dates}
    per_draw = len(dates) * len(ASSETS)
    rng = random.Random(seed)
    at_or_below_zero = 0
    for _ in range(replicates):
        picked = _draw_with_rng(dates, rng, block_dates)
        if sum(date_sums[date] for date in picked) / per_draw <= 0.0:
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
        "included_quarters": len(means),
        "positive_quarters": positive,
        "positive_fraction": positive / len(means) if means else None,
        "quarterly_pooled_means": means,
    }


def _rejection_key(exc: ValueError) -> str:
    message = str(exc)
    mapping = {
        "complete hourly grid missing": "incomplete-hourly-grid",
        "nonpositive price in required bar": "nonpositive-price-required-bar",
    }
    if message in mapping:
        return mapping[message]
    raise exc


def run_window(contract: dict, label: str, inventory: dict) -> dict:
    window = contract["windows"][label]
    minimum = contract["minimum_sample"][label]
    start_ms = parse_iso_ms(window["start_utc_inclusive"])
    end_ms = parse_iso_ms(window["end_utc_exclusive"])
    dates = expected_dates(start_ms, end_ms)

    complete: dict[str, dict[int, dict]] = {}
    rejected: dict[str, dict[str, int]] = {asset: {} for asset in ASSETS}
    for asset in ASSETS:
        bars_by_day, duplicated_days = load_asset_days(
            inventory, asset, start_ms, end_ms)
        observations: dict[int, dict] = {}
        for day_ms in dates:
            if day_ms in duplicated_days:
                key = "duplicate-timestamp-bar"
                rejected[asset][key] = rejected[asset].get(key, 0) + 1
                continue
            try:
                observations[day_ms] = daily_observation(
                    asset, day_ms, bars_by_day.get(day_ms, {}))
            except ValueError as exc:
                key = _rejection_key(exc)
                rejected[asset][key] = rejected[asset].get(key, 0) + 1
        complete[asset] = observations

    paired_dates = sorted(set(complete[ASSETS[0]]) & set(complete[ASSETS[1]]))
    rows = [complete[asset][day_ms] for day_ms in paired_dates for asset in ASSETS]
    per_asset_n = {
        asset: sum(row["asset"] == asset for row in rows) for asset in ASSETS
    }
    included_quarters = len({row["quarter"] for row in rows})
    sample = {
        "expected_dates": len(dates),
        "paired_dates": len(paired_dates),
        "pooled_observations": len(rows),
        "per_asset_observations": per_asset_n,
        "included_quarters": included_quarters,
        "rejected_asset_dates": rejected,
        "unpaired_complete_dates": {
            asset: len(set(complete[asset]) - set(paired_dates)) for asset in ASSETS
        },
        "minimum": minimum,
    }
    sample_ok = (
        len(paired_dates) >= minimum["paired_dates"]
        and len(rows) >= minimum["pooled"]
        and included_quarters >= minimum["included_quarters"]
        and all(per_asset_n[asset] >= minimum["per_asset"] for asset in ASSETS)
    )
    if not sample_ok:
        return {
            "window": window,
            "status": "inconclusive",
            "reason": "frozen minimum sample not met",
            "sample": sample,
            "gate_results": None,
        }

    pooled = {
        "mean_net": _mean([row["net"] for row in rows]),
        "mean_x": _mean([row["x"] for row in rows]),
    }
    per_asset = {}
    for asset in ASSETS:
        subset = [row for row in rows if row["asset"] == asset]
        per_asset[asset] = {
            "observations": len(subset),
            "mean_net": _mean([row["net"] for row in subset]),
            "mean_x": _mean([row["x"] for row in subset]),
        }
    stability = _quarter_stability(rows, "net")
    # Fresh identical RNG streams give both endpoints identical date draws.
    p_net = block_bootstrap_p(rows, "net")
    p_x = block_bootstrap_p(rows, "x")
    alpha = contract["loop_multiplicity"]["alpha_k"]

    gates = {
        "G1_net_magnitude": {
            "pooled_mean_net": pooled["mean_net"],
            "threshold": 0.0003,
            "pass": pooled["mean_net"] >= 0.0003,
        },
        "G2_net_asset_consistency": {
            "per_asset_mean_net": {
                asset: per_asset[asset]["mean_net"] for asset in ASSETS
            },
            "threshold": 0.0,
            "pass": all(per_asset[asset]["mean_net"] > 0.0 for asset in ASSETS),
        },
        "G3_quarter_stability": {
            "net": stability,
            "threshold": 0.60,
            # Exact rational comparison: positive/included >= 3/5.
            "pass": (stability["positive_quarters"] * 5
                     >= stability["included_quarters"] * 3),
        },
        "G4_concentration": {
            "pooled_mean_x": pooled["mean_x"],
            "per_asset_mean_x": {
                asset: per_asset[asset]["mean_x"] for asset in ASSETS
            },
            "threshold": 0.0,
            "pass": (pooled["mean_x"] > 0.0
                     and all(per_asset[asset]["mean_x"] > 0.0 for asset in ASSETS)),
        },
        "G5_dependence_robustness": {
            "one_sided_block_bootstrap_p_net": p_net,
            "one_sided_block_bootstrap_p_x": p_x,
            "alpha": alpha,
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": BOOTSTRAP_SEED,
            "block_dates": BOOTSTRAP_BLOCK_DATES,
            "pass": p_net < alpha and p_x < alpha,
        },
    }
    all_pass = all(gate["pass"] for gate in gates.values())
    return {
        "window": window,
        "status": "pass-on-this-sample" if all_pass else "falsified",
        "reason": None if all_pass else "at least one frozen gate is false",
        "sample": sample,
        "pooled": pooled,
        "per_asset": per_asset,
        "cost_model": {"taker_fee_per_side": TAKER_FEE,
                       "round_trip_fee_drag": 2 * TAKER_FEE,
                       "funding_crossed": False},
        "gate_results": gates,
        "gate_descriptions": contract["gates"],
    }


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
            args.primary_artifact.resolve().relative_to(WORKSPACE_ROOT))
        payload["primary_artifact_sha256"] = primary_sha
    write_json_new(args.output, payload)
    print(json.dumps({
        f"{args.window}_status": result["status"],
        "paired_dates": result["sample"]["paired_dates"],
        "pooled_observations": result["sample"]["pooled_observations"],
        "output": str(args.output),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
